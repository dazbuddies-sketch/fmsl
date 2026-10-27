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
import math, re

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


def prompt_for(item):
    v = item["variables"]
    return PROMPT.format(vars=", ".join(v), zeros=", ".join("0" for _ in v), f=item["f"])


def _num(s):
    s = s.strip().replace(" ", "")
    if s in ("+oo", "oo", "+inf", "inf", "+infinity", "infinity", "∞", "+∞"): return math.inf
    if s in ("-oo", "-inf", "-infinity", "-∞"): return -math.inf
    m = re.fullmatch(r"([+-]?\d+(?:\.\d+)?)/(\d+(?:\.\d+)?)", s)
    if m: return float(m.group(1)) / float(m.group(2))
    m = re.fullmatch(r"([+-]?)sqrt\((\d+(?:\.\d+)?)\)(?:/(\d+(?:\.\d+)?))?", s)
    if m:
        v = math.sqrt(float(m.group(2))) / (float(m.group(3)) if m.group(3) else 1.0)
        return -v if m.group(1) == "-" else v
    return float(s)


def parse_answer(text):
    """('LIMIT', value) | ('SET', [(lo, hi), ...]) | None"""
    clean = (text or "").replace("**", "").replace("__", "")
    lines = [l.strip() for l in clean.strip().splitlines() if "ANSWER:" in l.upper()]
    if not lines: return None
    line = lines[-1]
    m = re.search(r"ANSWER:\s*(LIMIT|SET)\s*(.+)$", line, re.IGNORECASE)
    if not m: return None
    kind, body = m.group(1).upper(), m.group(2).strip().rstrip(".").replace("`", "").replace("$", "")
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
    except (ValueError, ZeroDivisionError):
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
