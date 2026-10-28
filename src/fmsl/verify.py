"""Proven answers and answer checking, in plain data: the layer an AI assistant (or a grader,
or a test suite) calls instead of guessing.

    limit_at_origin("x^3*y/(x^6+y^2)")
        -> {"limit_exists": False, "set": "[-1/2, 1/2]", "proven": True, ...}
    check_claim("x^3*y/(x^6+y^2)", "LIMIT 0")
        -> {"correct": False, "missing": "[-1/2, 0) u (0, 1/2]", ...}

Expressions use the FMSL syntax (^ for powers; sin, cos, exp, abs, sign in one variable).
One variable: any expression FMSL's language handles. Two or three variables: rational functions.
Variables other than x, y, z are renamed internally; the point is always the origin.
"""
import ast, math, re
from fractions import Fraction as Fr

INF = math.inf
_NAMES = ("sin", "cos", "exp", "abs", "sign", "sqrt", "i", "pi")


# ------------------------------------------------------------------ helpers
def _variables(expr, variables=None):
    if variables:
        return list(variables)
    names = sorted(set(re.findall(r"[A-Za-z_]\w*", expr)) - set(_NAMES))
    if not names:
        raise ValueError("no variable found in the expression")
    return names


def _rename(expr, variables, target):
    out = expr
    tmp = {v: f"__v{k}__" for k, v in enumerate(variables)}
    for v, t in tmp.items(): out = re.sub(rf"\b{re.escape(v)}\b", t, out)
    for (v, t), new in zip(tmp.items(), target): out = out.replace(t, new)
    return out


def _f(v):
    """an engine value (Fraction, algebraic number, _Inf) as a float"""
    from fmsl.cluster import _Inf
    if isinstance(v, _Inf): return INF if v.sgn > 0 else -INF
    return float(v)


def _fmt(v):
    if v == INF: return "+oo"
    if v == -INF: return "-oo"
    fr = Fr(v).limit_denominator(1000)
    if abs(float(fr) - v) < 1e-12: return str(fr.numerator) if fr.denominator == 1 else f"{fr.numerator}/{fr.denominator}"
    return f"{v:.6g}"


def set_text(comps):
    return " u ".join("{" + _fmt(a) + "}" if a == b else f"[{_fmt(a)}, {_fmt(b)}]" for a, b in comps) or "(empty)"


# ------------------------------------------------------------------ the proven answer
def _one_variable(expr, var):
    from fmsl.compiler import run
    lines = []
    run(f"fluid f({var}) = {expr}\nat {var} -> 0: f\n", out=lines.append)
    txt = " ".join(lines)
    if "ERROR" in txt:
        raise ValueError(txt.split("ERROR", 1)[1].strip(" :"))
    proven = "UNPROVEN" not in txt and "not decided" not in txt
    conv = lambda s: {"+oo": INF, "-oo": -INF, "+inf": INF, "-inf": -INF}.get(s, None) if "oo" in s or "inf" in s else float(Fr(s))
    comps = None
    try:
        if "one-sided analysis" in txt:
            m = re.search(r"limit (\S+) on both sides", txt)
            if m: v = conv(m.group(1)); comps = [(v, v)]
            else:
                l = re.search(r"x->0-: (\S+?);", txt).group(1); r = re.search(r"x->0\+: (\S+?)(;| \[)", txt).group(1)
                comps = sorted({(conv(l), conv(l)), (conv(r), conv(r))})
        elif "ESSENTIAL" in txt:
            if "limit 0 on both sides" in txt: comps = [(0.0, 0.0)]
            elif "oscillates" in txt: comps = None
            else:
                side = lambda sd: conv(re.search(re.escape(sd) + r": ([^ ;(]+)", txt).group(1))
                comps = sorted({(side("x->0-"),) * 2, (side("x->0+"),) * 2})
        elif "removable" in txt:
            m = re.search(r"removable: value (\S+)", txt)
            v = 0.0 if "vanishes" in txt else float(Fr(m.group(1).rstrip(",")))
            comps = [(v, v)]
        elif "pole of order" in txt:
            m = re.search(r"pole of order (\d+).*\(lead ([^)]+)\)", txt)
            n, c = int(m.group(1)), float(Fr(m.group(2).strip()))
            right = INF if c > 0 else -INF
            left = right if n % 2 == 0 else -right
            comps = sorted({(left, left), (right, right)})
    except (AttributeError, ValueError, ZeroDivisionError):
        comps = None
    if comps is None: proven = False
    return comps, proven, txt.split(": f = ", 1)[-1].strip(), []


def _two_variables(expr):
    from fmsl.multivar import profile
    from fmsl.scaling import deciding_paths
    r = profile(expr)
    cl = r["cluster"]
    proven = bool(cl.get("exact")) and bool(cl.get("complete", True))
    comps = [(_f(a), _f(b)) for a, b in cl["signed_components"]]
    verdict, reasons = r["verdict"]
    witnesses = []
    try:
        for v, paths in deciding_paths(r):
            witnesses.append({"value": _fmt(_f(v)), "approached_along": paths[:3]})
    except Exception:                                   # witnesses are a convenience, never needed for the answer
        witnesses = []
    return comps, proven, verdict + ("; " + "; ".join(reasons) if reasons else ""), witnesses


def _three_variables(expr):
    from fmsl.space import profile3
    r = profile3(expr)
    proven = bool(r.get("exact", True)) and bool(r.get("complete"))
    comps = [(_f(a), _f(b)) for a, b in r["signed_components"]]
    return comps, proven, str(r.get("verdict", "")), []


def limit_at_origin(expr, variables=None):
    """The complete set of values expr approaches as the variables tend to 0, proven.

    Returns a dict: variables, limit_exists, limit (text or None), set (text), components (list of
    [lo, hi] floats on the extended line, inf for infinity), proven (bool), explanation, and, in two
    variables, witnesses: paths along which each end value is approached."""
    vs = _variables(expr, variables)
    if len(vs) == 1:
        comps, proven, expl, wit = _one_variable(_rename(expr, vs, ["x"]), "x")
    elif len(vs) == 2:
        comps, proven, expl, wit = _two_variables(_rename(expr, vs, ["x", "y"]))
    elif len(vs) == 3:
        comps, proven, expl, wit = _three_variables(_rename(expr, vs, ["x", "y", "z"]))
    else:
        raise ValueError("FMSL handles one, two or three variables")
    if len(vs) > 1:                                     # report paths in the caller's variable names
        back = dict(zip(["x", "y", "z"], vs))
        ren = lambda s: re.sub(r"\b([xyz])\b", lambda m: back.get(m.group(1), m.group(1)), s)
        expl = ren(expl)
        wit = [{"value": w["value"], "approached_along": [ren(p) for p in w["approached_along"]]} for w in wit]
    exists = comps is not None and len(comps) == 1 and comps[0][0] == comps[0][1] and not math.isinf(comps[0][0])
    return {"expression": expr, "variables": vs, "point": "origin",
            "limit_exists": exists if comps is not None else None,
            "limit": _fmt(comps[0][0]) if exists else None,
            "set": set_text(comps) if comps is not None else None,
            "components": [[a, b] for a, b in comps] if comps is not None else None,
            "proven": proven, "explanation": expl, "witnesses": wit}


# ------------------------------------------------------------------ checking a claimed answer
def _arith(n):
    if isinstance(n, ast.Constant) and isinstance(n.value, (int, float)): return float(n.value)
    if isinstance(n, ast.UnaryOp) and isinstance(n.op, (ast.UAdd, ast.USub)):
        v = _arith(n.operand); return -v if isinstance(n.op, ast.USub) else v
    if isinstance(n, ast.BinOp) and isinstance(n.op, (ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow)):
        a, b = _arith(n.left), _arith(n.right)
        if isinstance(n.op, ast.Add): return a + b
        if isinstance(n.op, ast.Sub): return a - b
        if isinstance(n.op, ast.Mult): return a * b
        if isinstance(n.op, ast.Div): return a / b
        if abs(b) > 10: raise ValueError("exponent too large")
        return a ** b
    if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id == "sqrt" and len(n.args) == 1:
        return math.sqrt(_arith(n.args[0]))
    raise ValueError("not a number")


def _num(s):
    s = s.strip().replace(" ", "")
    if s in ("+oo", "oo", "+inf", "inf", "∞", "+∞", "infinity", "+infinity"): return INF
    if s in ("-oo", "-inf", "-∞", "-infinity"): return -INF
    s = re.sub(r"√\s*(\d+(?:\.\d+)?)", r"sqrt(\1)", s).replace("√", "sqrt").replace("^", "**")
    s = re.sub(r"(\d)\s*(sqrt|\()", r"\1*\2", s)
    return _arith(ast.parse(s, mode="eval").body)


def parse_claim(text):
    """'LIMIT 3', '3', 'SET [-1/2, 1/2]', '[-oo, -2] u [2, +oo]', '{0} u [1, 2]', 'DNE'... -> list of (lo, hi)"""
    t = text.strip().replace("**", "^").replace("`", "").replace("$", "")
    t = re.sub(r"^(ANSWER:\s*)", "", t, flags=re.I)
    t = re.split(r"\s*(?:≈|~|\\approx)", t)[0].strip().rstrip(".")
    m = re.match(r"(LIMIT|SET)\s*(.*)$", t, re.I)
    if m: t = m.group(2).strip()
    if not re.search(r"[\[\(\{]", t):
        v = _num(t); return [(v, v)]
    out = []
    for part in re.split(r"\s*(?:\bu\b|\bU\b|∪|\\cup)\s*", t):
        mm = re.fullmatch(r"[\[\(]\s*(.+?)\s*,\s*(.+?)\s*[\]\)]", part.strip())
        if mm: out.append((_num(mm.group(1)), _num(mm.group(2)))); continue
        mm = re.fullmatch(r"\{\s*(.+?)\s*\}", part.strip())
        if mm:
            for v in mm.group(1).split(","): w = _num(v); out.append((w, w))
            continue
        raise ValueError(f"cannot read {part!r}")
    return sorted(out)


def _close(a, b, tol):
    if math.isinf(a) or math.isinf(b): return a == b
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def _minus(A, B, tol):
    """A minus B for lists of closed intervals on the extended line, as (lo, hi, lo_open, hi_open)
    pieces; slivers narrower than tol are dropped"""
    out = []
    for lo, hi in A:
        pieces = [(lo, hi, False, False)]
        for blo, bhi in B:
            nxt = []
            for plo, phi, lop, hop in pieces:
                if plo == phi:                          # a point: removed if B covers it
                    covered = (blo <= plo <= bhi) or _close(plo, blo, tol) or _close(plo, bhi, tol)
                    if not covered: nxt.append((plo, phi, lop, hop))
                    continue
                if bhi < plo or blo > phi or (_close(bhi, plo, tol) and bhi <= plo) or (_close(blo, phi, tol) and blo >= phi):
                    nxt.append((plo, phi, lop, hop)); continue
                if blo > plo and not _close(blo, plo, tol): nxt.append((plo, blo, lop, True))
                if bhi < phi and not _close(bhi, phi, tol): nxt.append((bhi, phi, True, hop))
            pieces = nxt
        out += pieces
    return out


def _pieces_text(P):
    return " u ".join("{" + _fmt(a) + "}" if a == b else ("(" if lo else "[") + f"{_fmt(a)}, {_fmt(b)}" + (")" if ho else "]")
                      for a, b, lo, ho in P)


def check_claim(expr, claim, variables=None, tol=1e-3):
    """Check a claimed answer (the limit, or the set of all limit values) against FMSL's proven one.

    Returns: correct (bool or None if FMSL could not prove the answer), truth, claimed, missing
    (values that are approached but the claim leaves out), extra (values the claim includes that
    are never approached), witnesses, and a one-line verdict."""
    truth = limit_at_origin(expr, variables)
    C = parse_claim(claim)
    res = {"expression": expr, "claim": claim, "claimed_set": set_text(C), "truth": truth["set"],
           "limit_exists": truth["limit_exists"], "proven": truth["proven"], "witnesses": truth["witnesses"]}
    if truth["components"] is None or not truth["proven"]:
        res.update(correct=None, missing=None, extra=None,
                   verdict="FMSL could not prove the answer for this expression; nothing is claimed.")
        return res
    T = [tuple(c) for c in truth["components"]]
    missing, extra = _minus(T, C, tol), _minus(C, T, tol)
    res.update(correct=not missing and not extra,
               missing=_pieces_text(missing) if missing else None, extra=_pieces_text(extra) if extra else None)
    if res["correct"]:
        res["verdict"] = f"Correct: the set of limit values is {truth['set']} (proven)."
    else:
        bits = []
        if extra: bits.append(f"these claimed values are never approached: {_pieces_text(extra)}")
        if missing: bits.append(f"these values are approached but missing from the claim: {_pieces_text(missing)}")
        res["verdict"] = f"Incorrect: the proven set is {truth['set']}; " + "; ".join(bits) + "."
    return res


# ------------------------------------------------------------------ the rate
def rate_at_origin(expr, variables=None):
    """How fast the limit is approached: the best theta with |f - L| <= C r^theta near the origin
    (proven, and attained along the reported arc). One or two variables."""
    from fractions import Fraction
    from .rate import rate, rate_1d
    vs = _variables(expr, variables)
    if len(vs) == 1:
        r = rate_1d(_rename(expr, vs, ["x"]), "x")
        text = r["text"].replace("|x|", f"|{vs[0]}|").replace("x ->", f"{vs[0]} ->")
    elif len(vs) == 2:
        r = rate(_rename(expr, vs, ["x", "y"]))
        back = dict(zip(["x", "y"], vs))
        text = re.sub(r"\b([xy])\b", lambda m: back[m.group(1)], r["text"])
    else:
        raise ValueError("the rate is computed in one or two variables")
    th = r.get("exponent")
    L = r.get("limit")
    return {"expression": expr, "variables": vs,
            "limit": None if L is None else (_fmt(float(L)) if not isinstance(L, str) else L),
            "exponent": None if th is None else (str(th) if isinstance(th, Fraction) else str(th)),
            "proven": r.get("proven", False), "explanation": text}
