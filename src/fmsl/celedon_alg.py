"""
fmsl.celedon_alg -- Celedon points inside zooms along curves that are IRRATIONAL in t.

fmsl.celedon handles zoom roots c0(t) in Q(t).  Here c0 is a root u of a polynomial h(u, t) with
no root in Q(t) (for example u^2 - t, u^2 - 2, u^2 - (t + 1)).  The substitution

        x = sigma * s^q,   y = s^p * (u + w)                          (a = p/q, s > 0)

gives a family whose coefficients are polynomials in (t, u), to be read on the curve h(u, t) = 0.
Every candidate condition of fmsl.moving.candidate_times is a polynomial P(t, u) built from
resultants in c; here these resultants are computed exactly (evaluation at several values of u,
interpolation, and a check at one more value), and u is then eliminated:

        R(t) = Res_u(P(t, u), h(u, t)).

Every t where P vanishes at some root u of h(., t) is a root of R.  So the real roots of R contain
every change along the irrational branch (and possibly conjugate-branch times: extra candidates
are harmless, fmsl.moving keeps only real changes of type).

A coefficient that vanishes ON the branch (h divides it) must be removed from the support before
the Newton polygon is read.  That is decided with a gcd over Q(t); if h factors so that the
coefficient vanishes on one part only, h is split and each part is treated separately.

Scope (reported, never hidden):
  * one irrational zoom level; a further zoom inside an irrational branch is reported unsupported;
  * end relations of the branch (boundary and critical values) are computed and returned, so
    arcs merging with these ends are found by fmsl.celedon.merge_candidate_times.
"""
from fractions import Fraction as Fr
from math import comb

from .algebraic import ptrim, pmul, padd, squarefree
from .cluster import _real_roots
from .plane import bgcd, bdivexact, resultant_d
from .moving import _tie_weights, _cdeg, _interp


# ----------------------------------------------------------- polynomials in (t, u)[c]
# a "tuc" polynomial is a dict {(t power, u power, c power): Fraction}
# a "tu" polynomial is a dict {(t power, u power): Fraction}

def _degc(P): return max((j for (_, _, j) in P), default=-1)
def _degu(P): return max((l for (_, l, _) in P), default=0)
def _degu2(P): return max((l for (_, l) in P), default=0)


def _mul3(P, Q):
    out = {}
    for (a, b, c), v in P.items():
        for (d, e, f), w in Q.items():
            k = (a + d, b + e, c + f)
            out[k] = out.get(k, 0) + v * w
    return {k: v for k, v in out.items() if v != 0}


def _sub3(P, Q):
    out = dict(P)
    for k, v in Q.items(): out[k] = out.get(k, 0) - v
    return {k: v for k, v in out.items() if v != 0}


def _dc3(P): return {(k, l, j - 1): v * j for (k, l, j), v in P.items() if j > 0}


def _at_u(P, u0):
    """tuc -> bivariate {(t power, c power)} at u = u0"""
    out = {}
    for (k, l, j), v in P.items():
        out[(k, j)] = out.get((k, j), 0) + v * u0 ** l
    return {kk: vv for kk, vv in out.items() if vv != 0}


def _lc_c_at_u(P, u0):
    d = _degc(P)
    return ptrim(_poly_t({(k, l): v for (k, l, j), v in P.items() if j == d}, u0))


def _poly_t(P2, u0):
    """tu -> polynomial in t at u = u0"""
    L = []
    for (k, l), v in P2.items():
        while len(L) <= k: L.append(Fr(0))
        L[k] += v * u0 ** l
    return ptrim(L)


def res_c(P, Q):
    """Res_c(P, Q) for tuc polynomials, as a tu polynomial (exact: evaluation in u at enough
       points where the leading coefficients in c survive, interpolation, one extra check)"""
    dP, dQ = _degc(P), _degc(Q)
    if dP <= 0 or dQ <= 0: return None
    bound = dP * _degu(Q) + dQ * _degu(P)
    us, rows, u0 = [], [], Fr(0)
    while len(us) < bound + 2:
        if _lc_c_at_u(P, u0) and _lc_c_at_u(Q, u0):
            us.append(u0); rows.append(ptrim(resultant_d(_at_u(P, u0), _at_u(Q, u0))))
        u0 = -u0 if u0 > 0 else -u0 + 1
        if abs(u0) > 10 * (bound + 5): raise ArithmeticError("no good evaluation points")
    dt = max((len(r) for r in rows), default=0)
    out = {}
    for k in range(dt):
        cs = _interp(us[:bound + 1], [r[k] if k < len(r) else Fr(0) for r in rows[:bound + 1]])
        for l, c in enumerate(cs):
            if c: out[(k, l)] = c
    # check at the extra point
    if _poly_t(out, us[-1]) != ptrim(rows[-1]): raise ArithmeticError("resultant interpolation check failed")
    return out


def norm(P2, h):
    """Res_u(P2, h) as a polynomial in t (P2, h: tu polynomials, u the eliminated variable)"""
    if not P2: return []
    if _degu2(P2) == 0:
        # P2 does not involve u: its roots are the conditions themselves
        L = []
        for (k, _), v in P2.items():
            while len(L) <= k: L.append(Fr(0))
            L[k] += v
        return ptrim(L)
    return ptrim(resultant_d(P2, h))


def _zero_on(P2, h):
    """how P2 (tu) sits on the curve h = 0: 'all' (h divides it), 'none', or a proper factor of h
       on which it vanishes (to split h)"""
    if not P2: return "all"
    g = bgcd(P2, h) if _degu2(P2) > 0 else ({(0, 0): Fr(1)} if P2 else {})
    dg = max((l for (_, l) in g), default=0)
    if dg == 0: return "none"
    if dg == _degu2(h): return "all"
    return g


# -------------------------------------------------------------------- the branch
def substitute_alg(N, D, a, sigma):
    """N, D in x^i y^j t^k -> dicts {(s power, w power): tu polynomial} under
       x = sigma s^q, y = s^p (u + w)"""
    p, q = a.numerator, a.denominator
    def sub(F):
        out = {}
        for (i, j, k), v in F.items():
            si = q * i + p * j
            sg = sigma ** i
            for e in range(j + 1):                     # (u + w)^j = sum C(j,e) u^(j-e) w^e
                key = (si, e)
                C = out.setdefault(key, {})
                kk = (k, j - e)
                C[kk] = C.get(kk, 0) + v * sg * comb(j, e)
        return {m: {kk: vv for kk, vv in C.items() if vv != 0} for m, C in out.items()
                if any(vv != 0 for vv in C.values())}
    return sub(N), sub(D)


def _initial3(C, b, sigma):
    """initial form (tuc) of {(i, j): tu} under s = sigma r^q', w = c r^p' (weight b)"""
    if not C: return {}
    wmin = min(i + b * j for (i, j) in C)
    out = {}
    for (i, j), P2 in C.items():
        if i + b * j != wmin: continue
        sg = sigma ** i
        for (k, l), v in P2.items():
            out[(k, l, j)] = out.get((k, l, j), 0) + sg * v
    return {kk: vv for kk, vv in out.items() if vv != 0}


def branch_candidates(N, D, a, sigma, h, label):
    """candidate polynomials in t along the irrational zoom branch u : h(u, t) = 0, the branch's
       end relations E(t, v) (norms over u), and what could not be covered.
       h: tu polynomial (dict {(t power, u power)}), squarefree in u."""
    polys, Es, unsupported = [], [], []
    CN, CD = substitute_alg(N, D, a, sigma)
    lc = {}
    du = _degu2(h)
    for (k, l), v in h.items():
        if l == du: lc[k] = lc.get(k, 0) + v
    L = ptrim([lc.get(k, Fr(0)) for k in range(max(lc) + 1)]) if lc else []
    if len(L) > 1: polys.append((L, f"the irrational zoom root escapes ({label})"))
    # remove coefficients that vanish on the branch; split h if needed
    for C in (CN, CD):
        for m in list(C):
            z = _zero_on(C[m], h)
            if z == "all":
                del C[m]
            elif z != "none":
                g = z
                return _split(N, D, a, sigma, h, g, label)
    if not CN or not CD:
        return polys, Es, unsupported
    # coefficients vanishing at some time on the branch
    for m, P2 in list(CN.items()) + list(CD.items()):
        r = norm(P2, h)
        if len(r) > 1: polys.append((r, f"coefficient of s^{m[0]} w^{m[1]} vanishes on the branch ({label})"))
    for b in _tie_weights(set(CN) | set(CD)):
        for sg in (1, -1):
            fN, fD = _initial3(CN, b, sg), _initial3(CD, b, sg)
            if not fN or not fD: continue
            side = "s>0" if sg > 0 else "s<0"
            g = _sub3(_mul3(_dc3(fN), fD), _mul3(fN, _dc3(fD)))
            jobs = [(fN, _dc3(fN), "numerator form has a repeated root", True),
                    (fD, _dc3(fD), "denominator form has a repeated root", True),
                    (fN, fD, "numerator and denominator forms share a root", True),
                    (g, _dc3(g), "critical points of the ratio collide", False),
                    (g, fD, "a critical point meets a pole", False)]
            for P, Q, why, zoom_if_zero in jobs:
                if not P or not Q or _degc(P) <= 0 or _degc(Q) <= 0: continue
                Pz, Qz = _drop_c3(P), _drop_c3(Q)
                if _degc(Pz) <= 0 or _degc(Qz) <= 0: continue
                r = res_c(Pz, Qz)
                if r is None: continue
                if _zero_on(r, h) != "none":
                    if zoom_if_zero:
                        unsupported.append(f"a further zoom inside the irrational branch {label} "
                                           f"(weight {b}, {side}): sampling only there")
                    continue
                R = norm(r, h)
                if len(R) > 1: polys.append((R, f"{why} (weight {b}, {side}) [on {label}]"))
            for P, nm in ((fN, "numerator"), (fD, "denominator")):
                for jj in (_degc(P), min(j for (_, _, j) in P)):
                    P2 = {(k, l): v for (k, l, j), v in P.items() if j == jj}
                    R = norm(P2, h)
                    if len(R) > 1: polys.append((R, f"{nm} form degenerates (weight {b}, {side}) [on {label}]"))
            Es.extend(_end_relations(fN, fD, g, h))
    return polys, Es, unsupported


def _drop_c3(P):
    m = min((j for (_, _, j) in P), default=0)
    return {(k, l, j - m): v for (k, l, j), v in P.items()} if m else P


def _split(N, D, a, sigma, h, g, label):
    h2 = bdivexact(h, g)
    out_p, out_E, out_u = [], [], []
    for part in (g, h2):
        p_, E_, u_ = branch_candidates(N, D, a, sigma, part, label)
        out_p += p_; out_E += E_; out_u += u_
    return out_p, out_E, out_u


def _end_relations(fN, fD, g, h):
    """end relations E(t, v) along the branch: boundary values (ratio of extreme coefficients)
       and critical values, with u eliminated (norms over u)"""
    out = []
    # boundary values: v * D_jj - N_jj = 0, for the extreme powers of c
    for jj in {_degc(fN), _degc(fD), min(j for (_, _, j) in fN), min(j for (_, _, j) in fD)}:
        Dj = {(k, l): v for (k, l, j), v in fD.items() if j == jj}
        Nj = {(k, l): v for (k, l, j), v in fN.items() if j == jj}
        if not Dj: continue
        E = _norm_linear_in_v(Dj, Nj, h)
        if E: out.append(E)
    # critical values: Res_c(g, v fD - fN) = 0, then Res_u(., h)
    gz = _drop_c3(g) if g else {}
    if gz and _degc(gz) > 0:
        try:
            E = _norm_critical(gz, fN, fD, h)
            if E: out.append(E)
        except ArithmeticError:
            pass
    return out


def _vpoints(n):
    return [Fr(k, 1) + Fr(1, 3) for k in range(n)]


def _norm_linear_in_v(Dj, Nj, h):
    """Res_u(v*Dj - Nj, h) as a dict {(t power, v power)} (v-degree <= deg_u h)"""
    m = _degu2(h)
    vs = _vpoints(m + 2)
    rows = []
    for v0 in vs:
        P2 = {}
        for kk, c in Dj.items(): P2[kk] = P2.get(kk, 0) + v0 * c
        for kk, c in Nj.items(): P2[kk] = P2.get(kk, 0) - c
        P2 = {kk: c for kk, c in P2.items() if c != 0}
        rows.append(norm(P2, h))
    return _interp_v(vs, rows, m)


def _norm_critical(g, fN, fD, h):
    """Res_u(Res_c(g, v fD - fN), h) as a dict {(t power, v power)}"""
    m = _degu2(h)
    dv = max(1, m) * _degc(g)
    vs = _vpoints(dv + 2)
    rows = []
    for v0 in vs:
        Q = _sub3({kk: v0 * c for kk, c in fD.items()}, fN)
        if _degc(Q) <= 0: return None
        r = res_c(g, Q)
        rows.append(norm(r, h) if r else [])
    return _interp_v(vs, rows, dv)


def _interp_v(vs, rows, deg):
    """rows[i] = polynomial in t at v = vs[i]; interpolate in v (degree <= deg), check the rest"""
    dt = max((len(r) for r in rows), default=0)
    E = {}
    for k in range(dt):
        cs = _interp(vs[:deg + 1], [r[k] if k < len(r) else Fr(0) for r in rows[:deg + 1]])
        for e, c in enumerate(cs):
            if c: E[(k, e)] = c
    for v0, r in zip(vs[deg + 1:], rows[deg + 1:]):
        val = []
        for (k, e), c in E.items():
            while len(val) <= k: val.append(Fr(0))
            val[k] += c * v0 ** e
        if ptrim(val) != ptrim(r): return None
    if not any(e >= 1 for (_, e) in E): return None
    return E
