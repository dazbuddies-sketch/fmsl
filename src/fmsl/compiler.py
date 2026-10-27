"""
fmsl.compiler -- the FMSL language: parser, analyzer, code generator.

Statements (one per line, '#' comments):
  context NAME = EXPR         fixed parameter (after Celedon's Meta-Value operator); exact constant
  fluid NAME(VAR) = EXPR      define a fluid value; transition point VAR -> 0
  at VAR -> 0: NAME           report the proven behavior of NAME at the point
  check EXPR                  three-valued zero test: proven_zero | proven_nonzero | unproven
  if NAME == 0 ...            REJECTED by the type checker (equality on fluid values is unsafe)
  emit python NAME            generate safe Python (series branch near 0 + direct formula)
  emit c NAME                 generate safe C, compile it with the system C compiler, validate it

Expressions: + - * / ^(integer)  sin cos exp  numbers (exact decimals)  i (imaginary unit)
             defined NAMEs, NAME(EXPR). Formulas containing i are compiled as complex-valued.
Engine: fmsl.engine (certified lazy Laurent series at 0). Numbers are exact rationals; pi/e/logs
are not allowed as constants because coefficients must stay rational.
"""
import cmath, math, re, sys, os, shutil, subprocess, tempfile
from fractions import Fraction as Fr
from . import engine as fms4
from .engine import V, Exp, Sin, Cos, Essential, EssFactor, infinity_flow
from .gauss import G
from .siletz import proven_series_error

# =============================== parser ===============================
TOK = re.compile(r"\s*(?:(\d+\.\d*|\d*\.\d+|\d+)|([A-Za-z_]\w*)|(\*\*|[-+*/^(),]))")

def tokenize(s):
    out, i = [], 0
    s = s.strip()
    while i < len(s):
        m = TOK.match(s, i)
        if not m or m.end() == i: raise SyntaxError(f"bad character at: {s[i:]!r}")
        num, name, op = m.groups()
        out.append(("num", num) if num else ("name", name) if name else ("op", "^" if op == "**" else op))
        i = m.end()
    return out

class Parser:
    def __init__(s, text): s.t, s.i = tokenize(text), 0
    def peek(s): return s.t[s.i] if s.i < len(s.t) else (None, None)
    def take(s, v=None):
        tok = s.peek()
        if v is not None and tok[1] != v: raise SyntaxError(f"expected {v!r}, got {tok[1]!r}")
        s.i += 1; return tok
    def parse(s):
        e = s.expr()
        if s.i != len(s.t): raise SyntaxError(f"unexpected {s.peek()[1]!r}")
        return e
    def expr(s):
        e = s.term()
        while s.peek()[1] in ("+", "-"):
            op = s.take()[1]; e = (op, e, s.term())
        return e
    def term(s):
        e = s.unary()
        while s.peek()[1] in ("*", "/"):
            op = s.take()[1]; e = (op, e, s.unary())
        return e
    def unary(s):
        if s.peek()[1] == "-": s.take(); return ("neg", s.unary())
        if s.peek()[1] == "+": s.take(); return s.unary()
        return s.power()
    def power(s):
        b = s.atom()
        if s.peek()[1] == "^":
            s.take()
            sign = -1 if s.peek()[1] == "-" else 1
            if s.peek()[1] in ("-", "+"): s.take()
            k, v = s.take()
            if k != "num" or "." in v: raise SyntaxError("exponent must be an integer")
            return ("pow", b, sign * int(v))
        return b
    def atom(s):
        k, v = s.take()
        if k == "num": return ("num", Fr(v))
        if k == "name":
            if v == "i" and s.peek()[1] != "(":
                return ("num", G(0, 1))
            if s.peek()[1] == "(":
                s.take(); a = s.expr(); s.take(")"); return ("call", v, a)
            return ("var", v)
        if v == "(":
            e = s.expr(); s.take(")"); return e
        raise SyntaxError(f"unexpected {v!r}")

FUNCS = ("sin", "cos", "exp", "abs", "sign")
SIDED = ("abs", "sign")        # functions analysed separately on each side of 0

def free_vars(e, defs):
    t = e[0]
    if t == "num": return set()
    if t == "var": return set() if e[1] in defs else {e[1]}
    if t == "call": return free_vars(e[2], defs)
    if t in ("neg",): return free_vars(e[1], defs)
    if t == "pow": return free_vars(e[1], defs)
    return free_vars(e[1], defs) | free_vars(e[2], defs)

def inline(e, defs, var):
    """Replace defined names: NAME -> body(var), NAME(arg) -> body(arg)."""
    t = e[0]
    if t == "num": return e
    if t == "var":
        if e[1] in defs:
            p, body = defs[e[1]]; return subst(body, p, ("var", var))
        return e
    if t == "call":
        a = inline(e[2], defs, var)
        if e[1] in defs:
            p, body = defs[e[1]]; return subst(body, p, a)
        if e[1] not in FUNCS: raise NameError(f"unknown function {e[1]!r}")
        return ("call", e[1], a)
    if t == "neg": return ("neg", inline(e[1], defs, var))
    if t == "pow": return ("pow", inline(e[1], defs, var), e[2])
    return (t, inline(e[1], defs, var), inline(e[2], defs, var))

def subst(e, p, a):
    t = e[0]
    if t == "num": return e
    if t == "var": return a if e[1] == p else e
    if t == "call": return ("call", e[1], subst(e[2], p, a))
    if t == "neg": return ("neg", subst(e[1], p, a))
    if t == "pow": return ("pow", subst(e[1], p, a), e[2])
    return (t, subst(e[1], p, a), subst(e[2], p, a))

# ============================ FMS semantics ============================
class NeedsSides(Exception):
    """the expression contains abs/sign: it must be analysed on each side of 0 separately"""


def _side_sign(f, side, fn):
    """the sign of f on the side `side` (+1: x > 0, -1: x < 0) of 0, exactly.
       Near 0, f = c x^k (1 + o(1)) with c != 0, so its sign there is sgn(c) * side^k."""
    if isinstance(f, EssFactor):
        raise Essential(f._undecided(f"inside {fn}(...)"))
    vv = f.val()
    if vv[0] == "zero":
        if vv[1]: return 0
        raise ValueError(f"{fn}(...) of a quantity indistinguishable from 0 through order {vv[2]}: "
                         f"its sign is not decided")
    k, c = vv[1], G(f.coeff(vv[1]))
    if not c.is_real(): raise ValueError(f"{fn}(...) of a complex quantity")
    sg = 1 if c.re > 0 else -1
    return sg * (side ** (k % 2))


def build(e, var, side=None):
    """AST -> certified fluid value (fms4) at var -> 0.  With abs/sign, `side` (+1 or -1) says
       which side of 0 is analysed; without a side, NeedsSides is raised."""
    t = e[0]
    if t == "call" and e[1] in SIDED:
        if side is None: raise NeedsSides(f"{e[1]}(...) is analysed on each side of 0: use 'at x -> 0'")
        f = build(e[2], var, side)
        sg = _side_sign(f, side, e[1])
        if e[1] == "sign": return V.rat([sg])
        return f * sg if sg != 0 else V.rat([0])
    if t == "num": return V.rat([e[1]])
    if t == "var":
        if e[1] != var: raise NameError(f"unknown name {e[1]!r}")
        return fms4.x
    if t == "neg": return -build(e[1], var, side)
    if t == "pow":
        b = build(e[1], var, side)
        return b ** e[2] if e[2] >= 0 else 1 / (b ** -e[2])
    if t == "call": return {"sin": Sin, "cos": Cos, "exp": Exp}[e[1]](build(e[2], var, side), keep=True)
    a, b = build(e[1], var, side), build(e[2], var, side)
    import operator
    return {"+": operator.add, "-": operator.sub, "*": operator.mul, "/": operator.truediv}[t](a, b)

def is_complex(e):
    t = e[0]
    if t == "num": return isinstance(e[1], G) and not e[1].is_real()
    if t == "var": return False
    if t in ("neg", "pow"): return is_complex(e[1])
    if t == "call": return is_complex(e[2])
    return is_complex(e[1]) or is_complex(e[2])

def _num_py(v):
    v = G(v) if not isinstance(v, G) else v
    if not v.is_real(): return repr(complex(v))
    r = v.re
    return repr(float(r)) if r.denominator != 1 else f"{r.numerator}.0"

def _num_c(v):
    v = G(v) if not isinstance(v, G) else v
    if not v.is_real(): return f"({float(v.re)!r} + {float(v.im)!r}*I)"
    r = v.re
    return repr(float(r)) if r.denominator != 1 else f"{r.numerator}.0"

def to_py(e, var, cx=False):
    t = e[0]
    if t == "num": return _num_py(e[1])
    if t == "var": return var
    if t == "neg": return f"(-{to_py(e[1], var, cx)})"
    if t == "pow": return f"({to_py(e[1], var, cx)})**{e[2]}"
    if t == "call": return f"{'cmath' if cx else 'math'}.{e[1]}({to_py(e[2], var, cx)})"
    return f"({to_py(e[1], var, cx)} {t} {to_py(e[2], var, cx)})"

def to_c(e, var, cx=False):
    t = e[0]
    if t == "num": return _num_c(e[1])
    if t == "var": return var
    if t == "neg": return f"(-{to_c(e[1], var, cx)})"
    if t == "pow":
        b, n = to_c(e[1], var, cx), e[2]
        if n == 0: return "1.0"
        pw = "cpow" if cx else "pow"
        prod = "*".join([f"({b})"] * abs(n)) if abs(n) <= 4 else f"{pw}({b}, {abs(n)})"
        return f"({prod})" if n > 0 else f"(1.0/({prod}))"
    if t == "call": return f"{'c' if cx else ''}{e[1]}({to_c(e[2], var, cx)})"
    return f"({to_c(e[1], var, cx)} {t} {to_c(e[2], var, cx)})"

def analyze(v):
    """Return (kind, text, proven) for a fluid value."""
    if isinstance(v, _SidedValue):
        return _analyze_sided(v)
    if isinstance(v, EssFactor):
        lim = v.limits()
        if lim is None:
            return ("essential", v._text(), True)
        return ("essential", "ESSENTIAL singularity: " + lim[2], True)
    try:
        vv = v.val()
    except Essential as ex:
        return ("essential", str(ex), True)
    if vv[0] == "zero":
        return ("zero", "C_e: identically zero [proven]" if vv[1]
                else f"zero through order {vv[2]} [UNPROVEN]", vv[1])
    k, c = vv[1], v.coeff(vv[1])
    proof = "[proven]" if vv[2] else "[exact leading term; zero-test class unproven]"
    if k > 0:  return ("removable", f"removable: value 0, vanishes to order {k} (lead {c}) {proof}", True)
    if k == 0: return ("removable", f"removable: value {c} {proof}", True)
    return ("pole", f"pole of order {-k}: {infinity_flow(k, c)} (lead {c}) {proof}", True)

# ============================ code generation ============================
NREF = 140   # terms in the high-accuracy reference series

def series_data(v):
    vv = v.val()
    lo = vv[1]
    coeffs = [v.coeff(k) for k in range(lo, lo + NREF)]
    # radius estimate (root test on exact coefficients)
    R = math.inf
    for i in range(20, NREF):
        k = lo + i; c = coeffs[i]
        if c != 0 and k > 0:
            ac = G(c).abs_upper()
            lg = (math.log2(ac.numerator) - math.log2(ac.denominator)) / k   # log2 |c|^(1/k)
            R = min(R, 2 ** (-lg)) if lg > 0 else R
    return lo, coeffs, R

def _tonum(v):
    v = G(v)
    return float(v.re) if v.is_real() else complex(v)

def ref_eval(lo, coeffs, xv):
    X = Fr(xv); acc = G(0)
    for c in reversed(coeffs): acc = acc * X + c
    return _tonum(acc * X ** lo) if lo >= 0 else _tonum(acc / X ** (-lo))

def series_eval(lo, coeffs, K, xv):
    acc = 0.0
    for c in reversed(coeffs[:K]): acc = acc * xv + _tonum(c)
    return acc * xv ** lo

def grid(limit):
    pts = []
    p = 0.0
    while p <= 12:
        mag = 10 ** (-p)
        if mag <= limit:
            pts += [mag, -mag, mag * 0.37, -mag * 0.73]
        p += 0.25
    return pts

def rel(a, b):
    if b == 0: return abs(a)
    if cmath.isnan(a) or cmath.isinf(a): return math.inf
    return abs(a - b) / abs(b)

def design(v, ast, var):
    """Pick switch point h and series length K that minimize measured max error."""
    lo, coeffs, R = series_data(v)
    limit = min(1.0, R / 2)
    cx = is_complex(ast)
    direct = eval(f"lambda {var}: {to_py(ast, var, cx)}", {"math": math, "cmath": cmath})
    pts = grid(limit)
    refs = {p: ref_eval(lo, coeffs, p) for p in pts}
    def derr(p):
        try: return rel(direct(p), refs[p])
        except (ZeroDivisionError, OverflowError, ValueError): return math.inf
    naive = max(derr(p) for p in pts)
    best = None
    for h in [limit, limit / 2, 0.25, 0.1, 0.05, 0.02, 0.01, 1e-3, 1e-4]:
        if h > limit: continue
        K = None
        for k in range(1, NREF - 20):
            if all(rel(series_eval(lo, coeffs, k, s), refs.get(s) or ref_eval(lo, coeffs, s)) < 2.3e-16
                   for s in (h, -h, h * .9, -h * .9)):  # ~1 ulp
                K = k; break
        if K is None: continue
        worst = max((rel(series_eval(lo, coeffs, K, p), refs[p]) if abs(p) < h else derr(p))
                    for p in pts)
        if best is None or (worst, K) < (best[0], best[2]): best = (worst, h, K)
    if best is None: raise ValueError("no series length reached full precision; formula too ill-conditioned near 0")
    worst, h, K = best
    # Siletz proof for the series zone; allow a few extra terms if that tightens the proof
    proven, K_used = None, K
    for k2 in list(range(K, K + 9)) + list(range(K + 10, K + 41, 4)):
        pb, _ = proven_series_error(v, lo, coeffs, k2, h, cx, abs(lo))
        if pb is not None and (proven is None or pb < proven / 2):
            proven, K_used = pb, k2
    if K_used != K:
        K = K_used
        worst = max((rel(series_eval(lo, coeffs, K, p), refs[p]) if abs(p) < h else derr(p)) for p in pts)
    return dict(cx=cx, lo=lo, coeffs=coeffs, R=R, limit=limit, naive=naive, worst=worst,
                h=h, K=K, npts=len(pts), proven=proven)

def proof_text(d, var):
    if d["proven"] is None:
        return f"Siletz proof: not available for this formula (series zone measured only)"
    return (f"Siletz proof: |rel. error| <= {float(d['proven']):.1e} for all 0 < |{var}| < {d['h']:g}"
            f" (series zone; includes truncation and floating-point rounding)")

def horner_py(lo, coeffs, K, var):
    s = _num_py(coeffs[K - 1])
    for c in reversed(coeffs[:K - 1]): s = f"{_num_py(c)} + {var}*({s})"
    if lo == 0: return s
    pw = "*".join([var] * abs(lo))
    return f"({pw})*({s})" if lo > 0 else f"({s})/({pw})"

def at_zero(kind, lo, c0):
    c0 = G(c0)
    if lo > 0: return "0j" if kind == "cx" else "0.0"
    if lo == 0: return f"complex({_num_py(c0)})" if kind == "cx" else _num_py(c0)
    if not c0.is_real(): return "complex(math.nan, math.nan)"   # complex pole: no single value
    if (-lo) % 2 == 0: return "math.inf" if c0 > 0 else "-math.inf"
    return "math.nan"

def emit_python(name, var, v, ast, d, text):
    lo, co, K, h = d["lo"], d["coeffs"], d["K"], d["h"]
    hdr = "import cmath, math" if d["cx"] else "import math"
    return (f"{hdr}\n\ndef {name}({var}):\n"
            f"    # FMSL: {text}\n"
            f"    # series branch |{var}| < {h:g}, {K} exact terms; measured max rel. error {d['worst']:.1e}\n"
            f"    # {proof_text(d, var)}\n"
            f"    if {var} == 0: return {at_zero('cx' if d['cx'] else None, lo, co[0])}\n"
            f"    if abs({var}) < {h!r}:\n"
            f"        return {horner_py(lo, co, K, var)}\n"
            f"    return {to_py(ast, var, d['cx'])}\n")

def emit_c(name, var, v, ast, d, text):
    lo, co, K, h = d["lo"], d["coeffs"], d["K"], d["h"]
    cx = d["cx"]
    s = _num_c(co[K - 1])
    for c in reversed(co[:K - 1]): s = f"{_num_c(c)} + {var}*({s})"
    pw = "*".join([var] * abs(lo))
    if lo > 0: s = f"({pw})*({s})"
    if lo < 0: s = f"({s})/({pw})"
    z = at_zero(None, lo, co[0])
    if z.startswith("complex("): z = "(NAN + NAN*I)"
    elif z.startswith("(") and "j" in z: z = _num_c(co[0])
    z = z.replace("math.inf", "INFINITY").replace("math.nan", "NAN")
    ty = "double complex" if cx else "double"
    inc = "#include <complex.h>\n#include <math.h>\n" if cx else "#include <math.h>\n"
    return (f"{inc}/* FMSL: {text}\n   series branch |{var}| < {h:g}, {K} exact terms\n"
            f"   {proof_text(d, var)} */\n"
            f"{ty} {name}(double {var}) {{\n"
            f"    if ({var} == 0.0) return {z};\n"
            f"    if (fabs({var}) < {h!r}) return {s};\n"
            f"    return {to_c(ast, var, cx)};\n}}\n")

def validate_c(name, var, csrc, d):
    """Compile the emitted C with the system compiler and measure it against the exact reference."""
    pts = grid(d["limit"])
    cc = shutil.which("gcc") or shutil.which("cc") or shutil.which("clang")
    if cc is None:
        return "no C compiler found; emitted code not validated"
    with tempfile.TemporaryDirectory() as tmp:
        src = os.path.join(tmp, "t.c"); exe = os.path.join(tmp, "t")
        if d["cx"]:
            main = csrc + "#include <stdio.h>\n#include <stdlib.h>\nint main(int c, char** v){" \
                   "for(int i=1;i<c;i++){double complex z=" + name + "(strtod(v[i],0));" \
                   "printf(\"%.17g %.17g\\n\", creal(z), cimag(z));}}\n"
        else:
            main = csrc + "#include <stdio.h>\n#include <stdlib.h>\nint main(int c, char** v){" \
                   "for(int i=1;i<c;i++) printf(\"%.17g\\n\", " + name + "(strtod(v[i],0)));}\n"
        open(src, "w").write(main)
        r = subprocess.run([cc, "-O2", "-o", exe, src, "-lm"], capture_output=True, text=True)
        if r.returncode: return f"C compile FAILED: {r.stderr[:300]}"
        out = subprocess.run([exe] + [repr(p) for p in pts], capture_output=True, text=True).stdout.split()
    if d["cx"]:
        vals = [complex(float(a), float(b)) for a, b in zip(out[0::2], out[1::2])]
    else:
        vals = [float(o) for o in out]
    worst = max(rel(o, ref_eval(d["lo"], d["coeffs"], p)) for o, p in zip(vals, pts))
    return f"compiled with {os.path.basename(cc)} -O2; measured max rel. error over {len(pts)} points: {worst:.1e}"

# ================================ runner ================================
def run(program, out=print):
    defs, values = {}, {}
    for lineno, raw in enumerate(program.splitlines(), 1):
        line = raw.split("#", 1)[0].strip()
        if not line: continue
        try:
            m = re.match(r"context\s+(\w+)\s*=\s*(.+)$", line)
            if m:
                name, body = m.groups()
                if name in ("i",) + FUNCS: raise NameError(f"{name!r} is reserved")
                ast = inline(Parser(body).parse(), defs, "_")
                fv = free_vars(ast, defs)
                if fv: raise NameError(f"context value must be a constant; unknown name(s): {', '.join(sorted(fv))}")
                if _has_call(ast):
                    raise ValueError("context value must be exact (no sin/cos/exp): coefficients must stay rational")
                defs[name] = (None, ast)
                out(f"[{lineno}] context {name} = {body.strip()}")
                continue
            m = re.match(r"fluid\s+(\w+)\s*\(\s*(\w+)\s*\)\s*=\s*(.+)$", line)
            if m:
                name, var, body = m.groups()
                ast = Parser(body).parse()
                fv = free_vars(ast, defs) - {var}
                if fv: raise NameError(f"unknown name(s): {', '.join(sorted(fv))}")
                full = inline(ast, defs, var)
                defs[name] = (var, full)
                values[name] = (var, full, None)
                out(f"[{lineno}] fluid {name}({var}) defined")
                continue
            m = re.match(r"at\s+(\w+)\s*->\s*0\s*:\s*(\w+)$", line)
            if m:
                var, name = m.groups()
                v = _value(values, name)
                out(f"[{lineno}] at {var} -> 0: {name} = {analyze(v)[1]}")
                continue
            m = re.match(r"check\s+(.+)$", line)
            if m:
                ast = Parser(m.group(1)).parse()
                fv = free_vars(ast, defs)
                var = fv.pop() if fv else "t"
                if fv: raise NameError("check expression must use a single variable")
                v = build(inline(ast, defs, var), var)
                vv = v.val()
                verdict = ("proven_zero" if vv[0] == "zero" and vv[1] else
                           f"unproven (zero through order {vv[2]}; no certificate)" if vv[0] == "zero" else
                           f"proven_nonzero (order {vv[1]})")
                out(f"[{lineno}] check {m.group(1).strip()}  ->  {verdict}")
                continue
            m = re.match(r"if\s+(.+?)\s*==\s*(.+)$", line)
            if m:
                out(f"[{lineno}] TYPE ERROR: '==' on fluid values is not allowed (it may be unprovable).\n"
                    f"        use: check {m.group(1)} - ({m.group(2).split()[0]})   "
                    f"and handle proven_zero / proven_nonzero / unproven")
                continue
            m = re.match(r"emit\s+(python|c)\s+(\w+)$", line)
            if m:
                lang, name = m.groups()
                var, ast, _ = values[name]
                v = _value(values, name)
                kind, text, proven = analyze(v)
                if kind == "sided":
                    out(f"[{lineno}] emit refused for {name}: code generation for abs/sign is not supported yet"); continue
                if kind in ("essential", "zero"):
                    out(f"[{lineno}] emit refused for {name}: {text}"); continue
                d = design(v, ast, var)
                if lang == "python":
                    src = emit_python(name, var, v, ast, d, text)
                    ns = {}; exec(src, ns)
                    out(f"[{lineno}] emit python {name}:\n" + _indent(src))
                    out(f"        accuracy (max rel. error vs exact, {d['npts']} pts, |{var}| <= {d['limit']:.3g}): "
                        f"naive formula {d['naive']:.1e}   FMSL {d['worst']:.1e}")
                    out(f"        {proof_text(d, var)}")
                else:
                    src = emit_c(name, var, v, ast, d, text)
                    out(f"[{lineno}] emit c {name}:\n" + _indent(src))
                    out(f"        {validate_c(name, var, src, d)}   (naive formula: {d['naive']:.1e})")
                    out(f"        {proof_text(d, var)}")
                continue
            raise SyntaxError("unknown statement")
        except Exception as ex:
            out(f"[{lineno}] ERROR: {type(ex).__name__}: {ex}")

def _has_call(e):
    t = e[0]
    if t in ("num", "var"): return False
    if t == "call": return True
    if t in ("neg", "pow"): return _has_call(e[1])
    return _has_call(e[1]) or _has_call(e[2])

def _value(values, name):
    var, ast, v = values[name]
    if v is None:
        try: v = build(ast, var)
        except Essential as ex: v = _EssentialValue(str(ex))
        except NeedsSides: v = _SidedValue(ast, var)
        values[name] = (var, ast, v)
    return v


class _SidedValue:
    """an expression with abs/sign: one fluid value per side of 0"""
    def __init__(s, ast, var):
        s.sides = {}
        for sd in (-1, 1):
            try: s.sides[sd] = build(ast, var, sd)
            except Essential as ex: s.sides[sd] = _EssentialValue(str(ex))
    def val(s): raise ValueError("abs/sign: analysed on each side of 0 separately")


def _one_side(v, sd):
    """the limit of v as x -> 0 from side sd: a printable value, or None if not decided"""
    if isinstance(v, EssFactor):
        lim = v.limits()
        return None if lim is None else str(lim[0] if sd < 0 else lim[1])
    try:
        vv = v.val()
    except Essential:
        return None
    if vv[0] == "zero": return "0" if vv[1] else None
    k, c = vv[1], G(v.coeff(vv[1]))
    if k > 0: return "0"
    if k == 0: return str(c)
    if not c.is_real(): return None
    sgn = (1 if c.re > 0 else -1) * (sd ** ((-k) % 2))
    return "+oo" if sgn > 0 else "-oo"


def _analyze_sided(v):
    l, r = _one_side(v.sides[-1], -1), _one_side(v.sides[1], 1)
    if l is None or r is None:
        return ("essential", f"one-sided analysis (abs/sign): x->0-: {l or 'not decided'}; "
                             f"x->0+: {r or 'not decided'} [UNPROVEN]", False)
    if l == r and "oo" not in l:
        return ("sided", f"limit {l} on both sides (one-sided analysis: abs/sign) [proven]", True)
    if l == r:
        return ("sided", f"one-sided analysis (abs/sign): x->0-: {l}; x->0+: {r} [proven]", True)
    return ("sided", f"one-sided analysis (abs/sign): x->0-: {l}; x->0+: {r}; no limit [proven]", True)

class _EssentialValue:
    def __init__(s, msg): s.msg = msg
    def val(s): raise Essential(s.msg)

def _indent(s): return "\n".join("        " + l for l in s.rstrip().splitlines())

def compile_text(program):
    """Run an FMSL program and return its report as a list of lines."""
    lines = []
    run(program, out=lines.append)
    return lines
