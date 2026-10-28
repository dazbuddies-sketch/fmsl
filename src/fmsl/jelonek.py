"""The non-properness set (Jelonek set) of a polynomial map of the real plane (roadmap #6).

For F = (f, g): R^2 -> R^2, the Jelonek set S_F is the set of points (a, b) that F approaches
along some sequence of points going to infinity:

        S_F = { lim F(p_k) : |p_k| -> infinity, F(p_k) converges }.

F is proper exactly over the complement of S_F. It matters for the Jacobian conjecture: by
Jelonek's theorem, a real polynomial map of the plane whose Jacobian determinant never vanishes is
a bijection if and only if its non-properness set is empty (Pinchuk's map, which shows that the
real version of the conjecture fails, has a nonempty one). S_F is empty or a closed curve made of
images of polynomial curves t -> (A(t), B(t)).

FMSL computes the real S_F with the same exact Newton-Puiseux zooming it uses at 0/0 points,
now at infinity. Every way of going to infinity is an arc
        x = s^-1 (or -s^-1),   y = c1 s^e1 + c2 s^e2 + ... ,   s -> 0+
(or the same with x and y exchanged, for arcs along which x stays bounded). Substituting, f and g
become Laurent polynomials in s; F stays bounded along the arc only if every negative power of s
cancels. Along y = ... + c s^e (generic c), the order of f is a piecewise linear function of e;
where it is negative, only special c (real roots of the initial form) can cancel it, and FMSL
zooms there; where the orders of f and g are >= 0, the limits are polynomials in the free
coefficient c, and (A(c), B(c)), c real, is a piece of S_F.

    from fmsl.jelonek import jelonek_set
    jelonek_set("x", "x*y")["text"]     # the line a = 0: (0, t), t real

Results are exact (rationals and real algebraic numbers); the zoom depth is capped, and a capped
search is reported as incomplete, never silently.
"""
from fractions import Fraction as Fr
from math import comb

MAX_DEPTH = 10


def _poly(text):
    from .multivar import parse_rational
    N, D = parse_rational(text)
    if set(D) != {(0, 0)}: raise ValueError("polynomials only")
    d = D[(0, 0)]
    return {k: v / d for k, v in N.items()}


def _m(P, b):
    return min(i + b * j for (i, j) in P) if P else None


def _slopes(P):
    out, pts = set(), list(P)
    for (i1, j1) in pts:
        for (i2, j2) in pts:
            if j2 > j1:
                b = Fr(i1 - i2, j2 - j1)
                m = _m(P, b)
                if len({j for (i, j) in pts if i + b * j == m}) >= 2: out.add(b)
    return out


def _initial(P, b):
    m = _m(P, b)
    terms = {}
    for (i, j), v in P.items():
        if i + b * j == m: terms[j] = terms.get(j, 0) + v
    return [terms.get(j, 0) for j in range(max(terms) + 1)]


def _zoom(P, p, q, c):
    R = {}
    for (i, j), v in P.items():
        e = q * i + p * j
        pw = [1]
        for _ in range(j): pw.append(pw[-1] * c)
        for k in range(j + 1):
            coef = v * comb(j, k) * pw[j - k]
            if coef != 0: R[(e, k)] = R.get((e, k), 0) + coef
    return {k: v for k, v in R.items() if v != 0}


def _peval(p, c):
    out = 0
    for v in reversed(p): out = out * c + v
    return out


def _trim(p):
    p = list(p)
    while p and p[-1] == 0: p.pop()
    return p


def _real_roots_nonzero(p):
    from .cluster import _real_roots
    p = _trim(p)
    if len(p) <= 1: return []
    return [r for r in _real_roots(p) if r != 0]


class _Incomplete(Exception): pass


def _search(Pf, Pg, depth, top, arc, out):
    """Pf, Pg: f and g in (s, w), s -> 0+, Laurent in s; record pieces of S_F in out"""
    if depth > MAX_DEPTH: raise _Incomplete("zoom deeper than the cap")
    # candidate slopes: corners of the Newton polygons, and where a term of f or g has order 0
    # (there the limit can depend on the free coefficient c)
    zc = {Fr(-i, j) for P in (Pf, Pg) for (i, j) in P if j > 0}
    sl = sorted(_slopes(Pf) | _slopes(Pg) | zc)
    if not top: sl = [b for b in sl if b > 0]
    # sample slopes: every edge slope, and one inside each gap (constant limits live there)
    lo_end = Fr(0) if not top else (sl[0] - 1 if sl else Fr(0))
    samples = list(sl)
    pts = ([lo_end] if top else [Fr(0)]) + sl
    for u, v in zip(pts, pts[1:]): samples.append((u + v) / 2)
    samples.append((sl[-1] if sl else (Fr(0) if not top else Fr(0))) + 1)
    if top and not sl: samples.append(Fr(-1))
    for b in sorted(set(samples)):
        if not top and b <= 0: continue
        mf, mg = _m(Pf, b), _m(Pg, b)
        If, Ig = _initial(Pf, b), _initial(Pg, b)
        if mf >= 0 and mg >= 0:
            A = If if mf == 0 else [0]
            B = Ig if mg == 0 else [0]
            out.append({"A": _trim(A) or [0], "B": _trim(B) or [0], "arc": arc + [("c", b)]})
        neg = [I for I, m in ((If, mf), (Ig, mg)) if m < 0]
        zero = [I for I, m in ((If, mf), (Ig, mg)) if m == 0]
        if neg:
            roots = _real_roots_nonzero(neg[0])
            roots = [c for c in roots if all(_peval(I, c) == 0 for I in neg[1:])]
        else:
            roots = []
            for I in zero:
                for c in _real_roots_nonzero(I):
                    if not any(c == r for r in roots): roots.append(c)
        p, q = b.numerator, b.denominator
        for c in roots:
            _search(_zoom(Pf, p, q, c), _zoom(Pg, p, q, c), depth + 1, False, arc + [(c, b)], out)


def _arcs_to_infinity(f, g, out):
    for swap in (False, True):
        for sigma in (1, -1):
            def chart(P):
                if swap: P = {(j, i): v for (i, j), v in P.items()}
                return {(-i, j): v * sigma ** i for (i, j), v in P.items()}
            pieces = []
            _search(chart(f), chart(g), 0, True, [], pieces)
            for pc in pieces:
                pc["chart"] = ("y" if swap else "x", sigma)
            out += pieces


def _same(p, q):
    p, q = _trim(p), _trim(q)
    return len(p) == len(q) and all(a == b for a, b in zip(p, q))


def _ptext(p, var="t"):
    from .optimize import _fmt
    p = _trim(p) or [0]
    terms = []
    for k, v in enumerate(p):
        if v == 0: continue
        cv = _fmt(v)
        if k == 0: terms.append(cv)
        else:
            co = "" if cv == "1" else "-" if cv == "-1" else f"{cv}*"
            terms.append(f"{co}{var}" + (f"^{k}" if k > 1 else ""))
    return " + ".join(terms).replace("+ -", "- ") or "0"


def _neg_t(p):
    return [v if k % 2 == 0 else -v for k, v in enumerate(p)]


def _union(ivs):
    from functools import cmp_to_key
    from .cluster import cmp
    ivs = sorted(ivs, key=cmp_to_key(lambda u, v: cmp(u[0], v[0])))
    out = []
    for lo, hi in ivs:
        if out and cmp(lo, out[-1][1]) <= 0:
            if cmp(hi, out[-1][1]) > 0: out[-1] = (out[-1][0], hi)
        else: out.append((lo, hi))
    return out


def _iv_text(ivs):
    from .optimize import _fmt
    def br(v, left):
        inf = getattr(v, "sgn", None) is not None
        return ("(" if left else ")") if inf else ("[" if left else "]")
    return " u ".join("{" + _fmt(lo) + "}" if lo == hi else f"{br(lo, True)}{_fmt(lo)}, {_fmt(hi)}{br(hi, False)}"
                      for lo, hi in ivs)


def jelonek_set(f_text, g_text):
    """The real non-properness set of F = (f, g): exact pieces of lines a = const or b = const (with
       the exact range of the other coordinate), parametrized curves (A(t), B(t)), t real, and points."""
    from .plane import image_1d
    from .cluster import cmp
    from .optimize import _fmt
    f, g = _poly(f_text), _poly(g_text)
    raw = []
    try:
        _arcs_to_infinity(f, g, raw)
        complete = True
    except _Incomplete:
        complete = False
    lines, curves, points = [], [], []
    for pc in raw:
        A, B = _trim(pc["A"]) or [0], _trim(pc["B"]) or [0]
        if len(A) <= 1 and len(B) <= 1:
            pt = (A[0], B[0])
            if not any(cmp(pt[0], q[0]) == 0 and cmp(pt[1], q[1]) == 0 for q in points): points.append(pt)
        elif len(A) <= 1 or len(B) <= 1:                  # a piece of a line a = const or b = const
            fixed, const, other = ("a", A[0], B) if len(A) <= 1 else ("b", B[0], A)
            L = next((L for L in lines if L["fixed"] == fixed and cmp(L["value"], const) == 0), None)
            if L is None:
                L = {"fixed": fixed, "value": const, "ivs": [], "arcs": []}; lines.append(L)
            L["ivs"].extend(image_1d(other, [Fr(1)])); L["arcs"].append(pc)
        elif not any((_same(A, C["A"]) and _same(B, C["B"])) or (_same(_neg_t(A), C["A"]) and _same(_neg_t(B), C["B"]))
                     for C in curves):
            curves.append({"A": A, "B": B, "chart": pc["chart"], "arc": pc["arc"]})
    line_list = [{"fixed": L["fixed"], "value": L["value"], "ranges": _union(L["ivs"]), "arcs": L["arcs"]} for L in lines]
    def covered(pt):
        for L in line_list:
            v, w = (pt[0], pt[1]) if L["fixed"] == "a" else (pt[1], pt[0])
            if cmp(v, L["value"]) == 0 and any(cmp(lo, w) <= 0 <= cmp(hi, w) for lo, hi in L["ranges"]): return True
        for C in curves:
            from .cluster import _real_roots
            ra = _trim([C["A"][0] - pt[0]] + list(C["A"][1:]))
            if len(ra) > 1 and any(_peval(C["B"], t) == pt[1] for t in _real_roots(ra)): return True
        return False
    points = [p for p in points if not covered(p)]
    parts = []
    for L in line_list:
        other = "b" if L["fixed"] == "a" else "a"
        whole = len(L["ranges"]) == 1 and str(L["ranges"][0][0]) in ("-oo",) or \
                (len(L["ranges"]) == 1 and getattr(L["ranges"][0][0], "sgn", 0) < 0 and getattr(L["ranges"][0][1], "sgn", 0) > 0)
        parts.append(f"the line {L['fixed']} = {_fmt(L['value'])}" + ("" if whole else f", {other} in {_iv_text(L['ranges'])}"))
    parts += [f"the curve (a, b) = ({_ptext(C['A'])}, {_ptext(C['B'])}), t real" for C in curves]
    parts += [f"the point ({_fmt(p[0])}, {_fmt(p[1])})" for p in points]
    text = ("empty: F is proper (preimages of bounded sets are bounded)" if not parts
            else "; ".join(parts)) + ("" if complete else "  [INCOMPLETE: zoom cap reached]")
    return {"lines": line_list, "curves": curves, "points": points, "complete": complete,
            "empty": not parts and complete, "text": text}
