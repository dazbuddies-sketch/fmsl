"""
fmsl.scaling -- the scaling map (Dimensional Unit Scaling, DUS, from MJ Celedon's FMS).

A 0/0 point is probed at every SCALE: along y = c|x|^a, the variable y is given the "dimension"
a relative to x (x ~ s, y ~ s^a as s -> 0).  At a given scale every monomial x^i y^j has size
s^(i + a j), so the smallest exponent dominates.  The scaling map lists

  * the REGIMES: ranges of a in which one monomial of N and one of D dominate; there
    f ~ (c_N / c_D) x^(i_N - i_D) y^(j_N - j_D) has size s^delta(a), delta = w_N(a) - w_D(a), so
    f -> 0 (delta > 0), f -> infinity (delta < 0), or f is of order 1 (delta = 0);
  * the BALANCE SCALES: the finitely many a where two or more monomials tie (the corners of the
    Newton polygons); only there can the leading terms interact, and that is where the
    interesting values of a fluid state come from;
  * which scale DECIDES each end of the fluid state: the paths (from the exact analysis of
    fmsl.multivar) on which that end value is attained or approached.

This is the rigorous form of DUS used throughout FMSL: weighted (quasi-homogeneous) scaling,
the standard tool of Newton polygons and tropical geometry.  The regimes are exact
(rational arithmetic on exponents); the attribution of the ends uses the exact fluid state.
"""
from fractions import Fraction as Fr
from functools import cmp_to_key

from .multivar import parse_rational, profile_polys
from .cluster import cmp, image_of_rational_function, _Inf, POS, NEG, _fmt


def _weights(P, a):
    return {m: m[0] + a * m[1] for m in P}


def _dominant(P, a):
    w = _weights(P, a)
    lo = min(w.values())
    return lo, sorted(m for m in P if w[m] == lo)


def balance_scales(N, D):
    """the scales a > 0 where two monomials of N, or two of D, tie at the minimal weight"""
    out = set()
    for P in (N, D):
        S = list(P)
        for u in range(len(S)):
            for v in range(u + 1, len(S)):
                (i1, j1), (i2, j2) = S[u], S[v]
                if j1 != j2:
                    a = Fr(i1 - i2, j2 - j1)
                    if a > 0 and len(_dominant(P, a)[1]) > 1: out.add(a)
    return sorted(out)


def _mono(m, c):
    i, j = m
    parts = [p for p in (("x" if i == 1 else f"x^{i}") if i else "", ("y" if j == 1 else f"y^{j}") if j else "") if p]
    body = "*".join(parts) or "1"
    if c == 1: return body
    if c == -1: return "-" + body
    return f"{_fmt(c)}*{body}" if parts else _fmt(c)


def _form_text(P, ms):
    return " + ".join(_mono(m, P[m]) for m in ms).replace("+ -", "- ")


def regimes(N, D):
    """the scaling map as a list of dicts: scale range, dominant terms, delta, behaviour"""
    bs = balance_scales(N, D)
    pts = [Fr(0)] + bs
    out = []
    for k, a0 in enumerate(bs + [None]):
        lo = pts[k]
        hi = a0
        # an interior scale of the open range (lo, hi)
        a = (lo + hi) / 2 if hi is not None else (lo + 1 if lo else Fr(1))
        if hi is None and bs: a = bs[-1] + 1
        if hi is not None and k == 0: a = hi / 2
        wN, mN = _dominant(N, a); wD, mD = _dominant(D, a)
        delta = wN - wD
        beh = "f -> 0" if delta > 0 else ("|f| -> oo" if delta < 0 else "f of order 1")
        if delta == 0 and mN == mD:
            beh = f"f -> {_fmt(Fr(N[mN[0]]) / D[mD[0]])}"
        out.append({"lo": lo, "hi": hi, "kind": "regime", "N": _form_text(N, mN), "D": _form_text(D, mD),
                    "delta": delta, "behaviour": beh, "slope": (mN[0][1] - mD[0][1]) if len(mN) == 1 == len(mD) else None})
        if a0 is not None:
            wN, mN = _dominant(N, a0); wD, mD = _dominant(D, a0)
            delta = wN - wD
            beh = "f -> 0" if delta > 0 else ("|f| -> oo" if delta < 0 else "balance: f of order 1, depends on c")
            # where the denominator's form vanishes, a finer scale (a zoom) decides
            if len(mD) > 1:
                form = [Fr(0)] * (1 + max(m[1] for m in mD))
                for m in mD: form[m[1]] += D[m]
                from .cluster import _real_roots
                from .algebraic import ptrim
                rts = [r for r in _real_roots(ptrim(form)) if cmp(r, 0) != 0] if len(ptrim(form)) > 1 else []
                if rts:
                    beh += (f"; D's leading form vanishes at c = {', '.join(_vtext(r) for r in rts)}: "
                            f"a finer scale decides there")
            out.append({"a": a0, "kind": "balance", "N": _form_text(N, mN), "D": _form_text(D, mD),
                        "delta": delta, "behaviour": beh})
    return out


def _entry_set(e):
    """the closed set of values that a profile entry contributes, as (lo, hi) pairs"""
    k = e.get("kind")
    if k == "value": return [(e["value"], e["value"])]
    if k == "inf": return [(POS, POS) if s > 0 else (NEG, NEG) for s in (e.get("signs") or {1, -1})]
    if k == "all": return [(NEG, POS)]
    if k == "depends":
        n, d = e["R"]
        try:
            from .algebraic import pgcd, pdivmod, ptrim
            g = pgcd(n, d); n2, d2 = pdivmod(n, g)[0], pdivmod(d, g)[0]
            if len(ptrim(n2)) <= 1 and len(ptrim(d2)) <= 1:
                v = (ptrim(n2) or [Fr(0)])[0] / d2[0]; return [(v, v)]
            return image_of_rational_function(n2, d2)
        except (ArithmeticError, NotImplementedError, TypeError):
            return []
    return []


def deciding_paths(res):
    """for each end of the fluid state: the paths on which it is attained or approached"""
    cl = res["cluster"]
    comps = cl.get("signed_components") or []
    ends = []
    for lo, hi in comps:
        ends.append(lo)
        if cmp(lo, hi) != 0: ends.append(hi)
    out = []
    for v in ends:
        who = []
        for e in res["entries"]:
            for lo, hi in _entry_set(e):
                if cmp(lo, v) == 0 or cmp(hi, v) == 0:
                    who.append(e.get("path", "?")); break
        seen = []
        for w in who:
            if w not in seen: seen.append(w)
        # the most specific paths: zooms ('exactly', '+ c...') and balance scales come first
        def rank(p):
            return (0 if ("exactly" in p or " + " in p) else 1 if "c*|x|" in p else 2)
        seen.sort(key=rank)
        best = [p for p in seen if rank(p) == rank(seen[0])] if seen else []
        out.append((v, best))
    return out


def _vtext(v):
    if isinstance(v, _Inf): return "+oo" if v.sgn > 0 else "-oo"
    return _fmt(v)


def scaling_map(text):
    N, D = parse_rational(text)
    res = profile_polys(N, D)
    return {"regimes": regimes(N, D), "deciding": deciding_paths(res),
            "fluid_state": res["cluster"].get("signed_text")}


def describe_scaling(text):
    r = scaling_map(text)
    out = [f"f = {text}", "scaling map (y ~ c|x|^a; size of f ~ |x|^delta):"]
    for g in r["regimes"]:
        if g["kind"] == "regime":
            rng = f"{_fmt(g['lo'])} < a < {_fmt(g['hi'])}" if g["hi"] is not None else f"a > {_fmt(g['lo'])}"
            if g["lo"] == 0 and g["hi"] is not None: rng = f"0 < a < {_fmt(g['hi'])}"
            out.append(f"  {rng:16s} N ~ {g['N']},  D ~ {g['D']}:  {g['behaviour']}")
        else:
            out.append(f"  a = {_fmt(g['a']):12s} BALANCE  N ~ {g['N']},  D ~ {g['D']}:  {g['behaviour']}")
    out.append(f"fluid state: {r['fluid_state']}")
    for v, paths in r["deciding"]:
        out.append(f"  the end {_vtext(v)} is decided on: " + ("; ".join(paths[:3]) if paths else "(not attributed)")
                   + (f" (+{len(paths) - 3} more)" if len(paths) > 3 else ""))
    return "\n".join(out)
