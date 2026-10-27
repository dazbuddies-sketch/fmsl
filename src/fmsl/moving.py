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
            "sig": _signature(res),
            # a whole-line answer cannot be missing anything, even if some branch was unresolved
            "exact": cl.get("exact", False) and (cl.get("complete", True) or bool(cl.get("whole"))),
            "comps": cl.get("signed_components")}


def moving_state(N, D, t0):
    """all limits as (x, y, t) -> (0, 0, t0) together (three-variable engine, z = t - t0)"""
    N3, D3 = _shift_time(N, t0), _shift_time(D, t0)
    r = profile3_polys(N3, D3)
    return {"text": r["signed_text"], "verdict": r["verdict"], "complete": r["complete"]}


def _samples_between(lo, hi):
    """three rational times strictly inside (lo, hi); lo/hi None = unbounded"""
    if lo is None and hi is None: return [Fr(-1), Fr(0), Fr(1)]
    from math import floor, ceil
    if lo is None:
        h = floor(_bounds(hi)[0]) if isinstance(hi, A) else floor(hi)
        h = h if h < hi else h - 1
        return [Fr(h - 2), Fr(h - 1), Fr(h)]
    if hi is None:
        l = ceil(_bounds(lo)[1]) if isinstance(lo, A) else ceil(lo)
        l = l if l > lo else l + 1
        return [Fr(l), Fr(l + 1), Fr(l + 2)]
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
        ends = []
        if mid and mid.get("comps") is not None:
            try:
                ends = endpoint_formulas(N, D, pts[1], mid["comps"])
            except (ArithmeticError, NotImplementedError, TypeError):
                ends = []
        intervals.append({"lo": lo, "hi": hi, "sample": pts[1],
                          "snapshot": mid["text"] if mid else "undefined (D = 0)",
                          "exact": bool(mid and mid["exact"]), "ends": ends})
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


def _tpoly_text(L):
    """polynomial in t (ascending Fractions) as text"""
    terms = []
    for k in range(len(L) - 1, -1, -1):
        c = L[k]
        if not c: continue
        cs = str(int(c)) if c.denominator == 1 else f"({c})"
        mon = "" if k == 0 else ("t" if k == 1 else f"t^{k}")
        if mon:
            body = mon if abs(c) == 1 else f"{str(abs(int(c))) if c.denominator == 1 else f'({abs(c)})'}*{mon}"
        else:
            body = str(abs(int(c))) if c.denominator == 1 else f"{abs(c)}"
        terms.append(("-" if c < 0 else "+", body))
    if not terms: return "0"
    out = ("-" if terms[0][0] == "-" else "") + terms[0][1]
    for sg, body in terms[1:]: out += f" {sg} {body}"
    return out


def _solve_text(E):
    """E(t, v) = 0 as text; solved for v when E is linear in v"""
    if max(e for (_, e) in E) == 1:
        L1 = [Fr(0)] * (1 + max([k for (k, e) in E if e == 1] or [0]))
        L0 = [Fr(0)] * (1 + max([k for (k, e) in E if e == 0] or [0]))
        for (k, e), c in E.items(): (L1 if e == 1 else L0)[k] += c
        L1, L0 = ptrim(L1), ptrim([-c for c in L0])
        if len(L1) == 1:
            return "v = " + _tpoly_text([c / L1[0] for c in L0])
        return f"v = ({_tpoly_text(L0)}) / ({_tpoly_text(L1)})"
    return _E_text(E) + "  (v = the end)"


def describe_moving(text, path=None):
    r = moving_profile(text, path)
    out = [f"f = {text}" + (f",  target at ({path[0]}, {path[1]})" if path else "")]
    out.append("snapshots (fluid state at a fixed time):")
    for iv in r["intervals"]:
        lo = "-oo" if iv["lo"] is None else _tstr(iv["lo"])
        hi = "+oo" if iv["hi"] is None else _tstr(iv["hi"])
        more = "".join(f";  at t = {t}: {txt}" for t, txt in iv.get("more", [])[:2])
        out.append(f"  {lo} < t < {hi}:  at t = {iv['sample']}: {iv['snapshot']}{more}")
        forms = []
        for v, E in iv.get("ends", []):
            if E:
                txt = _solve_text(E)
                if txt not in forms: forms.append(txt)
        if forms:
            out.append(f"      ends move with t as:  " + ";  ".join(forms))
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


# ------------------------------------------------------------------ endpoint formulas
# A finite end of a snapshot is (at the top level of the analysis) either a critical value of
# the ratio R_t(c) = fN(t, c) / fD(t, c) of initial forms, or a limit of R_t at c -> 0 or oo.
# A critical value v satisfies E(t, v) = Res_c(R' numerator, v fD - fN) = 0.  E is computed
# exactly by evaluation at rational times and interpolation in t: its degree in t is at most
# the Sylvester bound B, so B + 1 exact values determine it; two further times are checked.

def _interp(ts, vals):
    """exact Lagrange interpolation: polynomial (ascending list) through (ts[k], vals[k])"""
    from .algebraic import padd, pmul
    out = []
    for k, (tk, yk) in enumerate(zip(ts, vals)):
        if not yk: continue
        num, den = [Fr(1)], Fr(1)
        for m, tm in enumerate(ts):
            if m == k: continue
            num = pmul(num, [-tm, Fr(1)]); den *= (tk - tm)
        out = padd(out, [c * yk / den for c in num])
    return ptrim(out)


def _endpoint_polys(N, D):
    """candidate polynomials E(t, v) as dicts {(t power, v power)}"""
    CN, CD = _coeff_polys(N), _coeff_polys(D)
    out = []
    for a in _tie_weights(set(CN) | set(CD)):
        for sigma in (1, -1):
            fN, fD = _initial_form(CN, a, sigma), _initial_form(CD, a, sigma)
            if not fN or not fD: continue
            for jj in {_cdeg(fN), _cdeg(fD), min(j for (_, j) in fN), min(j for (_, j) in fD)}:
                E = {}
                for (k, j), v in fD.items():
                    if j == jj: E[(k, 1)] = E.get((k, 1), 0) + v
                for (k, j), v in fN.items():
                    if j == jj: E[(k, 0)] = E.get((k, 0), 0) - v
                E = {kk: vv for kk, vv in E.items() if vv}
                if any(e == 1 for (_, e) in E): out.append(E)
            g = _bsub(_bmul(_dc(fN), fD), _bmul(fN, _dc(fD)))
            if not g or _cdeg(g) == 0: continue
            degt = lambda P: max((k for (k, _) in P), default=0)
            # bound on deg_t of Res_c(g, v fD - fN)
            dQ = max(_cdeg(fN), _cdeg(fD)); tQ = max(degt(fN), degt(fD))
            B = _cdeg(g) * tQ + dQ * degt(g)
            ts = [Fr(k, 1) + Fr(1, 7) for k in range(B + 3)]
            rows = []
            for t0 in ts:
                gt = {}
                for (k, j), v in g.items(): gt[(0, j)] = gt.get((0, j), 0) + v * t0 ** k
                Q = {}
                for (k, j), v in fD.items(): Q[(1, j)] = Q.get((1, j), 0) + v * t0 ** k
                for (k, j), v in fN.items(): Q[(0, j)] = Q.get((0, j), 0) - v * t0 ** k
                gt = {kk: vv for kk, vv in _merge(gt).items() if vv}
                Q = {kk: vv for kk, vv in Q.items() if vv}
                rows.append(ptrim(resultant_d(Q, gt)) if gt and Q else [])   # polynomial in v
            dv = max((len(r) for r in rows), default=0)
            E = {}
            for e in range(dv):
                coeffs_t = _interp(ts[:B + 1], [r[e] if e < len(r) else Fr(0) for r in rows[:B + 1]])
                for k, c in enumerate(coeffs_t):
                    if c: E[(k, e)] = c
            ok = all(_evalE(E, t0) == ptrim(r) for t0, r in zip(ts[B + 1:], rows[B + 1:]))
            if E and ok: out.append(E)
    return out


def _merge(P):
    out = {}
    for k, v in P.items(): out[k] = out.get(k, 0) + v
    return out


def _evalE(E, t0):
    """E(t0, v) as a polynomial in v"""
    dv = 1 + max((e for (_, e) in E), default=-1)
    L = [Fr(0)] * dv
    for (k, e), c in E.items(): L[e] += c * t0 ** k
    return ptrim(L)


def _primitive_E(E):
    """remove the content in t (a factor depending on t only) and scale to integers"""
    from math import gcd
    from .algebraic import pgcd
    byv = {}
    for (k, e), c in E.items():
        L = byv.setdefault(e, [])
        while len(L) <= k: L.append(Fr(0))
        L[k] += c
    cont = None
    for L in byv.values():
        cont = ptrim(L) if cont is None else pgcd(cont, ptrim(L))
    from .algebraic import pdivmod
    if cont and len(cont) > 1:
        byv = {e: pdivmod(ptrim(L), cont)[0] for e, L in byv.items()}
    E2 = {(k, e): c for e, L in byv.items() for k, c in enumerate(L) if c}
    den = 1
    for c in E2.values(): den = den * c.denominator // gcd(den, c.denominator)
    E2 = {k: c * den for k, c in E2.items()}
    g = 0
    for c in E2.values(): g = gcd(g, int(c))
    return {k: c / g for k, c in E2.items()} if g > 1 else E2


def _E_text(E):
    terms = []
    for (k, e), c in sorted(E.items(), key=lambda kv: (-kv[0][1], -kv[0][0])):
        c = int(c) if c.denominator == 1 else c
        mon = "*".join(p for p in (("t" if k == 1 else f"t^{k}") if k else "", ("v" if e == 1 else f"v^{e}") if e else "") if p)
        if mon:
            co = "" if c == 1 else ("-" if c == -1 else f"{c}*")
            terms.append(f"{co}{mon}")
        else:
            terms.append(str(c))
    s = " + ".join(terms).replace("+ -", "- ")
    return s + " = 0"


def endpoint_formulas(N, D, t0, comps):
    """for the finite ends of the snapshot at the rational time t0: an exact relation
       E(t, v) = 0 satisfied by that end (v) for t in this interval, when one is found"""
    from .algebraic import peval
    Es = [_primitive_E(E) for E in _endpoint_polys(N, D)]
    out = []
    ends = []
    for lo, hi in comps:
        for v in ((lo,) if cmp(lo, hi) == 0 else (lo, hi)):
            if not isinstance(v, _Inf): ends.append(v)
    for v in ends:
        hit = None
        for E in sorted(Es, key=lambda E: (max(e for (_, e) in E), len(E))):
            p = _evalE(E, t0)
            if len(p) > 1 and peval(p, v) == 0:
                hit = E; break
        out.append((v, hit))
    return out
