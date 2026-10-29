"""Grading for LimitBench.

A model must end its reply with one line:
    ANSWER: LIMIT <value>          (the limit exists and equals <value>)
    ANSWER: SET <set>              (no limit; <set> is the set of all limit values)
Endpoints are accepted within 1e-3 relative (so 4 significant digits suffice).
<set> is a union of closed intervals and points on the extended line, e.g. "[-1/2, 1/2]",
"[-oo, -3] u [3, +oo]", "{0} u [1, 2]"; numbers as integers, fractions or decimals; +oo, -oo.

Scores per item:
  verdict  -- right about whether the limit exists
  exact    -- the limit value, or the whole set, is right (endpoints within 1e-3 relative)
"""
import ast, math, re

PROMPT = """Let f({vars}) = {f}.
Determine the behaviour of f as ({vars}) -> ({zeros}): find the set of ALL limit values, that is,
every L (a real number, +oo or -oo) such that f(p_n) -> L for some sequence of points p_n -> ({zeros})
at which f is defined. If this set is a single real number, the limit exists.

Explain your reasoning, then end with exactly one final line in one of these two forms:
ANSWER: LIMIT <value>
ANSWER: SET <set>
where <set> is a union of closed intervals and points, for example [-1/2, 1/2] or [-oo, -3] u [3, +oo]
or {{0}} u [1, 2], using +oo and -oo for infinity, and numbers as integers, fractions, sqrt(...) or
decimals with at least 4 significant digits."""


PROMPT_INF = PROMPT.replace("as ({vars}) -> ({zeros})", "as |({vars})| -> oo (the point moves arbitrarily far from the origin, in any manner)") \
                   .replace("p_n -> ({zeros})", "p_n with |p_n| -> oo")


def prompt_for(item):
    v = item["variables"]
    if item.get("point") == "oo":
        return PROMPT_INF.format(vars=", ".join(v), zeros="", f=item["f"])
    text = PROMPT.format(vars=", ".join(v), zeros=", ".join("0" for _ in v), f=item["f"])
    side = item.get("side")
    if side:
        cond = " and ".join(f"{k} {'>' if s == '+' else '<'} 0" for k, s in side.items())
        text = text.replace("at which f is defined.", f"at which f is defined, restricted to {cond}.") \
                   .replace(f"as ({', '.join(v)}) -> ({', '.join('0' for _ in v)}):",
                            f"as ({', '.join(v)}) -> ({', '.join('0' for _ in v)}) with {cond}:")
    return text


def _num(s):
    s = s.strip().replace(" ", "")
    if s in ("+oo", "oo", "+inf", "inf", "+infinity", "infinity", "∞", "+∞"): return math.inf
    if s in ("-oo", "-inf", "-infinity", "-∞"): return -math.inf
    s = re.sub(r"√\s*(\d+(?:\.\d+)?)", r"sqrt(\1)", s).replace("√", "sqrt").replace("^", "**")
    s = re.sub(r"(\d)\s*(sqrt|\()", r"\1*\2", s)              # implicit product: 2sqrt(3), 2(1+x)
    return _arith(ast.parse(s, mode="eval").body)


def _arith(n):
    """Evaluate a closed-form number safely: + - * / ** and sqrt(...), nothing else."""
    if isinstance(n, ast.Constant) and isinstance(n.value, (int, float)): return float(n.value)
    if isinstance(n, ast.UnaryOp) and isinstance(n.op, (ast.UAdd, ast.USub)):
        v = _arith(n.operand); return -v if isinstance(n.op, ast.USub) else v
    if isinstance(n, ast.BinOp) and isinstance(n.op, (ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow)):
        a, b = _arith(n.left), _arith(n.right)
        if isinstance(n.op, ast.Pow) and abs(b) > 10: raise ValueError("exponent too large")
        return {ast.Add: a + b, ast.Sub: a - b, ast.Mult: a * b}.get(type(n.op)) if not isinstance(n.op, (ast.Div, ast.Pow)) \
            else (a / b if isinstance(n.op, ast.Div) else a ** b)
    if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id == "sqrt" and len(n.args) == 1:
        return math.sqrt(_arith(n.args[0]))
    raise ValueError("not a number")


def parse_answer(text):
    """('LIMIT', value) | ('SET', [(lo, hi), ...]) | None"""
    clean = (text or "").replace("**", "").replace("__", "")
    lines = [l.strip() for l in clean.strip().splitlines() if "ANSWER:" in l.upper()]
    if not lines: return None
    line = lines[-1]
    m = re.search(r"ANSWER:\s*(LIMIT|SET)\s*(.+)$", line, re.IGNORECASE)
    if not m: return None
    kind, body = m.group(1).upper(), m.group(2).strip().replace("`", "").replace("$", "")
    body = re.split(r"\s*(?:≈|~|\\approx)", body)[0].strip().rstrip(".")      # drop a trailing "≈ decimals"
    try:
        if kind == "LIMIT":
            return ("LIMIT", _num(body))
        comps = []
        for part in re.split(r"\s*(?:u|U|∪|\\cup)\s*", body):
            part = part.strip()
            mm = re.fullmatch(r"[\[\(]\s*(.+?)\s*,\s*(.+?)\s*[\]\)]", part)
            if mm: comps.append((_num(mm.group(1)), _num(mm.group(2)))); continue
            mm = re.fullmatch(r"\{\s*(.+?)\s*\}", part)
            if mm: v = _num(mm.group(1)); comps.append((v, v)); continue
            return None
        return ("SET", comps)
    except (ValueError, ZeroDivisionError, SyntaxError, TypeError, OverflowError):
        return None


def _close(a, b, tol=1e-3):
    if math.isinf(a) or math.isinf(b): return a == b
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def grade(item, answer_text):
    """{'parsed': bool, 'verdict': bool, 'exact': bool}"""
    a = parse_answer(answer_text)
    truth = [tuple(c) for c in item["fluid_state_float"]]
    truth_limit = item["limit_exists"]
    if a is None: return {"parsed": False, "verdict": False, "exact": False}
    if a[0] == "LIMIT":
        says_limit = not math.isinf(a[1])
        comps = [(a[1], a[1])]
    else:
        comps = sorted(a[1])
        says_limit = len(comps) == 1 and _close(comps[0][0], comps[0][1]) and not math.isinf(comps[0][0])
    verdict = says_limit == truth_limit
    exact = len(comps) == len(truth) and all(_close(x[0], y[0]) and _close(x[1], y[1])
                                              for x, y in zip(sorted(comps), sorted(truth)))
    return {"parsed": True, "verdict": verdict, "exact": exact}
