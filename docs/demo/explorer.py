"""Glue between the Fluid State Explorer page and FMSL (runs in the browser with Pyodide).

analyze(text) returns a dict for the page: the exact fluid state from FMSL, and, for comparison,
what the common "try a few paths" method concludes (limits along standard paths, computed
exactly here from the same rational function)."""
import re
from fractions import Fraction as Fr

from fmsl.multivar import profile, parse_rational
from fmsl.scaling import describe_scaling
from fmsl.space import profile3, parse3
from fmsl.moving import describe_moving


def _fmt(v):
    if v is None: return "no limit"
    if v == "+oo" or v == "-oo": return v
    return str(v) if v.denominator == 1 else f"{v.numerator}/{v.denominator}"


def _limit_along(terms_N, terms_D):
    """limit as s -> 0+ of (sum c s^e) / (sum c s^e), exact"""
    def lead(T):
        T = {e: c for e, c in T.items() if c != 0}
        if not T: return None
        e = min(T); return e, T[e]
    a, b = lead(terms_N), lead(terms_D)
    if b is None: return None                       # denominator identically 0 on the path
    if a is None: return Fr(0)
    (en, cn), (ed, cd) = a, b
    if en > ed: return Fr(0)
    if en == ed: return Fr(cn) / cd
    return "+oo" if cn / cd > 0 else "-oo"


def _path2(N, D, fx, fy):
    """substitute x = fx(s), y = fy(s), each (coefficient, power) or None for 0"""
    def sub(P):
        out = {}
        for (i, j), c in P.items():
            if (i and fx is None) or (j and fy is None): continue
            cx, ex = fx if i else (1, 0)
            cy, ey = fy if j else (1, 0)
            e = (ex * i if i else 0) + (ey * j if j else 0)
            out[e] = out.get(e, 0) + c * (Fr(cx) ** i) * (Fr(cy) ** j)
        return out
    return _limit_along(sub(N), sub(D))


PATHS2 = [("along y = 0", (1, 1), None), ("along x = 0", None, (1, 1)),
          ("along y = x", (1, 1), (1, 1)), ("along y = -x", (1, 1), (-1, 1)),
          ("along y = 2x", (1, 1), (2, 1)), ("along y = x^2", (1, 1), (1, 2)),
          ("along y = -x^2", (1, 1), (-1, 2)), ("along x = y^2", (1, 2), (1, 1))]


def path_test(text):
    N, D = parse_rational(text)
    rows = [(name, _path2(N, D, fx, fy)) for name, fx, fy in PATHS2]
    vals = {(_fmt(v)) for _, v in rows}
    concl = (f"all paths give {next(iter(vals))}: the path test concludes the limit is {next(iter(vals))}"
             if len(vals) == 1 else "the paths disagree: the path test concludes there is no limit")
    return [(n, _fmt(v)) for n, v in rows], concl


def analyze(text):
    text = text.strip()
    names = set(re.findall(r"\b[a-z]\b", text))
    if "sqrt" in text or "r" in names:
        from fmsl.radial import radial_limit
        r = radial_limit(text)
        return {"mode": "dist", "state": r["set"] or "not decided", "complete": bool(r["proven"]),
                "verdict": ("the limit exists: " + r["set"].strip("{}")) if r["limit_exists"] else
                           "no single limit: these are all the values approached"}
    if "t" in names:
        return {"mode": "moving", "report": describe_moving(text)}
    if "z" in names:
        r = profile3(text)
        return {"mode": "3d", "state": r["signed_text"], "verdict": r["verdict"],
                "complete": bool(r["complete"])}
    r = profile(text)
    cl = r["cluster"]
    v = r["verdict"]
    verdict = v[0] if isinstance(v, tuple) else v
    rows, concl = path_test(text)
    return {"mode": "2d", "state": cl.get("signed_text", cl.get("text")), "verdict": verdict,
            "exact": bool(cl.get("exact", False)), "complete": bool(cl.get("complete", True)),
            "paths": rows, "path_conclusion": concl, "scaling": describe_scaling(text)}
