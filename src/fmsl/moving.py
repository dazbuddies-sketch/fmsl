"""
fmsl.moving -- fluid states of a moving target.

A family f_t(x, y) = N(x, y, t) / D(x, y, t) (coefficients polynomial in the time t) has a
0/0 point at the origin for each t, or at a moving point p(t) = (p1(t), p2(t)) given as a path.
Two questions:

  * the SNAPSHOT at time t0: the fluid state of f_t0 at the origin (two-variable engine);
  * the MOVING state at t0: all limits of f(x, y, t) as (x, y, t) -> (0, 0, t0) together,
    i.e. the target is approached while it moves.  This is a three-variable problem: the
    three-variable engine (fmsl.space) with z = t - t0.

The snapshot changes character only at finitely many SPECIAL TIMES.  This module computes

  1. candidate special times, exactly: real roots of
       - every coefficient (a polynomial in t) of N and D  (the Newton polygons change there),
       - for every weight a where monomials tie: the discriminants and mutual resultants of the
         initial forms in c (y = c|x|^a), and of the critical-point polynomial of their ratio;
     all as exact polynomials in t (resultants over Q[t]);
  2. the snapshot at sample times in every interval between candidates, and a check that its
     type (the components and which ends are infinite or collapse) is the same at three points
     of each interval -- any change not explained by a candidate is located by bisection and
     reported as an approximate change point, never hidden;
  3. at every candidate where the type does change: the snapshot and the moving state.

What is proven and what is not: every reported fluid state is exact.  Candidate times are exact
algebraic numbers.  Completeness of the list of special times is proven for changes visible in
the top-level Newton data; changes that arise only deeper in the zoom are detected by sampling
(reported as "change near ..."), which is evidence, not proof.
"""
import re
from fractions import Fraction as Fr
from functools import cmp_to_key

from .space import parse3, profile3_polys, _cancel_monomial_content
from .multivar import profile_polys
from .plane import resultant_d
from .algebraic import A, ptrim, pderiv, squarefree
from .cluster import cmp, _real_roots, _rational_between, _bounds, set_unique, POS, NEG, _Inf


def parse_moving(text, path=None):
    """N, D as dicts {(i, j, k): coeff} for x^i y^j t^k.  path: (p1, p2) polynomials in t
       (strings): the 0/0 point moves to (p1(t), p2(t))."""
    s = text
    if path:
        p1, p2 = path
        s = re.sub(r"\bx\b", f"(x+({p1}))", s)
        s = re.sub(r"\by\b", f"(y+({p2}))", s)
    if re.search(r"\bz\b", s): raise ValueError("use x, y and the time t")
    s = re.sub(r"\bt\b", "z", s)
    return parse3(s)


def _at_time(P, t0):
    """P(x, y, t0) as a two-variable dict"""
    out = {}
    for (i, j, k), v in P.items():
        out[(i, j)] = out.get((i, j), 0) + v * (t0 ** k if k else 1)
    return {k: v for k, v in out.items() if v != 0}


def _shift_time(P, t0):
    """P(x, y, z + t0): the time measured from t0"""
    from math import comb
    out = {}
    for (i, j, k), v in P.items():
        for e in range(k + 1):
            c = v * comb(k, e) * (t0 ** (k - e) if k - e else 1)
            out[(i, j, e)] = out.get((i, j, e), 0) + c
    return {k: v for k, v in out.items() if v != 0}


def _coeff_polys(P):
    """{(i, j): coefficient polynomial in t (ascending list)}"""
    out = {}
    for (i, j, k), v in P.items():
        L = out.setdefault((i, j), [])
        while len(L) <= k: L.append(Fr(0))
        L[k] += v
    return {m: ptrim(L) for m, L in out.items() if ptrim(L)}


def _tie_weights(S):
    """weights a > 0 (y ~ |x|^a) where two monomials of the support S tie"""
    ws = set()
    S = list(S)
    for u in range(len(S)):
        for w in range(u + 1, len(S)):
            (i1, j1), (i2, j2) = S[u], S[w]
            if j1 != j2:
                a = Fr(i1 - i2, j2 - j1)
                if a > 0: ws.add(a)
    return sorted(ws)


def _initial_form(C, a, sigma):
    """minimal-weight part of sum C[(i,j)](t) x^i y^j under x = sigma s^q, y = c s^p, as a
       bivariate dict {(t power, c power)}"""
    if not C: return {}
    wmin = min(i + a * j for (i, j) in C)
    out = {}
    for (i, j), L in C.items():
        if i + a * j != wmin: continue
        sg = sigma ** i
        for k, v in enumerate(L):
            if v: out[(k, j)] = out.get((k, j), 0) + sg * v
    return {k: v for k, v in out.items() if v != 0}


def _dc(P): return {(k, j - 1): v * j for (k, j), v in P.items() if j > 0}

def _bmul(P, Q):
    out = {}
    for (a, b), u in P.items():
        for (c, d), w in Q.items():
            out[(a + c, b + d)] = out.get((a + c, b + d), 0) + u * w
    return {k: v for k, v in out.items() if v != 0}

def _bsub(P, Q):
    out = dict(P)
    for k, v in Q.items(): out[k] = out.get(k, 0) - v
    return {k: v for k, v in out.items() if v != 0}

def _cdeg(P): return max((j for (_, j) in P), default=0)


def candidate_times(N, D):
    """exact real candidate special times (sorted), with the reason each one appears"""
    CN, CD = _coeff_polys(N), _coeff_polys(D)
    polys = []                                           # (polynomial in t, reason)
    for (m, L) in list(CN.items()) + list(CD.items()):
        if len(L) > 1: polys.append((L, f"coefficient of x^{m[0]} y^{m[1]} vanishes"))
    for a in _tie_weights(set(CN) | set(CD)):
        for sigma in (1, -1):
            fN, fD = _initial_form(CN, a, sigma), _initial_form(CD, a, sigma)
            g = _bsub(_bmul(_dc(fN), fD), _bmul(fN, _dc(fD))) if fN and fD else {}
            side = "x>0" if sigma > 0 else "x<0"
            jobs = [(fN, _dc(fN), "numerator form has a repeated root"),
                    (fD, _dc(fD), "denominator form has a repeated root"),
                    (fN, fD, "numerator and denominator forms share a root"),
                    (g, _dc(g), "critical points of the ratio collide"),
                    (g, fD, "a critical point meets a pole")]
            for P, Q, why in jobs:
                if P and Q and _cdeg(P) > 0 and _cdeg(Q) > 0:
                    r = ptrim(resultant_d(P, Q))
                    if len(r) > 1: polys.append((r, f"{why} (weight {a}, {side})"))
            # leading / trailing coefficients in c (roots escape to c = 0 or infinity)
            for P, nm in ((fN, "numerator"), (fD, "denominator")):
                if not P: continue
                for jj in (_cdeg(P), min(j for (_, j) in P)):
                    L = [Fr(0)] * (1 + max(k for (k, j) in P if j == jj))
                    for (k, j), v in P.items():
                        if j == jj: L[k] += v
                    if len(ptrim(L)) > 1: polys.append((ptrim(L), f"{nm} form degenerates (weight {a}, {side})"))
    found = []                                           # (time, reasons)
    for L, why in polys:
        sq = squarefree(L)
        if len(sq) <= 1: continue
        for r in _real_roots(sq):
            for item in found:
                if cmp(item[0], r) == 0: item[1].add(why); break
            else:
                found.append([r, {why}])
    found.sort(key=cmp_to_key(lambda u, v: cmp(u[0], v[0])))
    return [(r, sorted(w)) for r, w in found]


def _kind(v):
    if isinstance(v, _Inf): return "+oo" if v.sgn > 0 else "-oo"
    return "fin"


def _signature(res):
    """the TYPE of a fluid state: verdict, and for each component which ends are infinite and
       whether it is a single point"""
    cl = res["cluster"]
    comps = cl.get("signed_components") or []
    v = res["verdict"]
    verdict = v[0] if isinstance(v, tuple) else v
    verdict = verdict.split(" (")[0].split(" =")[0]
    sig = tuple((_kind(lo), _kind(hi), cmp(lo, hi) == 0) for lo, hi in comps)
    return (verdict, sig)


def snapshot(N, D, t0):
    """the fluid state of f_t0 at the origin (two-variable engine); None if D(., ., t0) = 0"""
    N2, D2 = _at_time(N, t0), _at_time(D, t0)
    if not D2: return None
    if not N2:
        return {"text": "{0}", "verdict": "LIMIT", "sig": ("LIMIT", (("fin", "fin", True),)), "exact": True}
    res = profile_polys(N2, D2)
    cl = res["cluster"]
    return {"text": cl.get("signed_text", cl.get("text")), "verdict": res["verdict"],
            "sig": _signature(res), "exact": cl.get("exact", False) and cl.get("complete", True)}


def moving_state(N, D, t0):
    """all limits as (x, y, t) -> (0, 0, t0) together (three-variable engine, z = t - t0)"""
    N3, D3 = _shift_time(N, t0), _shift_time(D, t0)
    r = profile3_polys(N3, D3)
    return {"text": r["signed_text"], "verdict": r["verdict"], "complete": r["complete"]}


def _samples_between(lo, hi):
    """three rational times strictly inside (lo, hi); lo/hi None = unbounded"""
    if lo is None and hi is None: return [Fr(-1), Fr(0), Fr(1)]
    if lo is None:
        h = _bounds(hi)[0] if isinstance(hi, A) else hi
        return [h - 3, h - 2, h - 1] if not isinstance(hi, A) else [h - 3, h - 2, h - 1]
    if hi is None:
        l = _bounds(lo)[1] if isinstance(lo, A) else lo
        return [l + 1, l + 2, l + 3]
    m = _rational_between(lo, hi)
    m1 = _rational_between(lo, m)
    m2 = _rational_between(m, hi)
    return [m1, m, m2]


def _bisect_change(N, D, a, b, sa, sb, tol=Fr(1, 10 ** 6)):
    """locate a type change between rational times a < b (types sa != sb)"""
    while b - a > tol:
        m = (a + b) / 2
        s = snapshot(N, D, m)
        sm = s["sig"] if s else None
        if sm == sa: a = m
        else: b, sb = m, sm
    return a, b


def moving_profile(text, path=None, joint=True):
    """the full moving-target analysis of a family f(x, y, t).  Returns a dict with
       'intervals' (lo, hi, sample time, snapshot text), 'special' (time, reasons, snapshot,
       moving state), 'unexplained' (approximate change points found only by sampling)."""
    N, D = parse_moving(text, path)
    cands = candidate_times(N, D)
    times = [c[0] for c in cands]
    bounds = [None] + times + [None]
    intervals, unexplained = [], []
    sig_of_interval = []
    for k in range(len(bounds) - 1):
        lo, hi = bounds[k], bounds[k + 1]
        pts = _samples_between(lo, hi)
        snaps = [snapshot(N, D, p) for p in pts]
        sigs = [s["sig"] if s else None for s in snaps]
        for u in range(len(pts) - 1):
            if sigs[u] != sigs[u + 1]:
                a, b = _bisect_change(N, D, pts[u], pts[u + 1], sigs[u], sigs[u + 1])
                unexplained.append((a, b))
        mid = snaps[1]
        intervals.append({"lo": lo, "hi": hi, "sample": pts[1],
                          "snapshot": mid["text"] if mid else "undefined (D = 0)",
                          "exact": bool(mid and mid["exact"])})
        sig_of_interval.append(sigs[1])
    special = []
    for k, (t0, why) in enumerate(cands):
        try:
            s = snapshot(N, D, t0)
        except (ArithmeticError, NotImplementedError, TypeError) as exc:
            s = {"text": f"not computed ({type(exc).__name__})", "sig": "?", "exact": False}
        left, right = sig_of_interval[k], sig_of_interval[k + 1]
        here = s["sig"] if s else None
        if left == right == here:
            continue                                       # nothing changes at this candidate
        entry = {"t": t0, "why": why, "snapshot": s["text"] if s else "undefined (D = 0)"}
        if joint:
            try:
                m = moving_state(N, D, t0)
                entry["moving"] = m["text"] + ("" if m["complete"] else " (partial)")
            except (ArithmeticError, NotImplementedError, TypeError) as exc:
                entry["moving"] = f"not computed ({type(exc).__name__})"
        special.append(entry)
    # merge neighbouring intervals of the same type that no special time separates
    merged = []
    sp_times = [e["t"] for e in special]
    for k, iv in enumerate(intervals):
        if merged and sig_of_interval[k] == merged[-1][1] and iv["lo"] is not None \
                and not any(cmp(t, iv["lo"]) == 0 for t in sp_times) \
                and not any(a < (_bounds(iv["lo"])[1] if isinstance(iv["lo"], A) else iv["lo"]) < b
                            for a, b in unexplained):
            merged[-1][0]["hi"] = iv["hi"]
            merged[-1][0].setdefault("more", []).append((iv["sample"], iv["snapshot"]))
        else:
            merged.append([dict(iv), sig_of_interval[k]])
    return {"intervals": [m[0] for m in merged], "special": special, "unexplained": unexplained,
            "candidates": len(cands)}


def _tstr(t):
    if t is None: return "oo"
    return f"{float(t):.6g} (algebraic)" if isinstance(t, A) else str(t)


def describe_moving(text, path=None):
    r = moving_profile(text, path)
    out = [f"f = {text}" + (f",  target at ({path[0]}, {path[1]})" if path else "")]
    out.append("snapshots (fluid state at a fixed time):")
    for iv in r["intervals"]:
        lo = "-oo" if iv["lo"] is None else _tstr(iv["lo"])
        hi = "+oo" if iv["hi"] is None else _tstr(iv["hi"])
        more = "".join(f";  at t = {t}: {txt}" for t, txt in iv.get("more", [])[:2])
        out.append(f"  {lo} < t < {hi}:  at t = {iv['sample']}: {iv['snapshot']}{more}")
    if r["special"]:
        out.append("special times (the snapshot changes character):")
        for s in r["special"]:
            out.append(f"  t = {_tstr(s['t'])}:  snapshot {s['snapshot']};  moving target {s.get('moving', '-')}")
            kinds = sorted({re.sub(r" \(weight.*\)$", "", w) for w in s["why"]})
            out.append(f"      why: {'; '.join(kinds[:3])}" + (f" (+{len(kinds) - 3} more)" if len(kinds) > 3 else ""))
    for a, b in r["unexplained"]:
        out.append(f"  change near t in ({float(a):.6g}, {float(b):.6g}) not explained by the top-level "
                   f"candidates (found by sampling)")
    return "\n".join(out)
