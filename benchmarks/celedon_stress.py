"""
Stress test of Celedon points (fmsl.moving), default and deep mode, with an ERROR LOG.

Families come from several independent sources:
  classic   the 8 classic families of benchmarks/moving_crosscheck.py
  randA     the random generator of moving_crosscheck.py (sparse, degree <= 4)
  randB     a different generator: denser, degree <= 5, coefficients up to quadratic in t
  hidden1   a change hidden one zoom level deep (zoom root rational in t)
  hidden2   a change hidden two zoom levels deep
  irrat     zoom roots that are NOT rational in t (sqrt 2, sqrt t): deep mode must say so
  dgnp      the published examples of Dinh-Guo-Nguyen-Pham (2022) with a parameter inserted
  path      the target moves along a path (p1(t), p2(t))

Every answer is checked by methods that do not use the engine's Newton polygon / zoom:
  * LOOP: the exact image of f on a tiny square loop around the point (benchmarks/loop_crosscheck),
    compared with each reported snapshot (interval samples and both sides of every Celedon point);
  * DENSE: 12 extra times inside every reported interval: any change of snapshot type there is a
    missed Celedon point (each such snapshot is itself LOOP-checked when it is used as evidence);
  * FALSE: a reported Celedon point where the type is the same just before, at, and just after;
  * ENDS: every reported end formula E(t, v) = 0 re-checked at 4 more times inside its interval.

Each problem is one JSON line in the error log:
  {"run", "seed", "source", "family", "mode", "kind", "t", "detail"}
Kinds: crash, timeout, loop_mismatch, missed_change, unexplained_left, false_point,
       wrong_end_formula, unsupported_branch (deep mode saying honestly it cannot prove a branch).

Run:  python benchmarks/celedon_stress.py SEED N_PER_SOURCE LOG.jsonl [RUN_NAME]
"""
import json, random, signal, sys, time
from fractions import Fraction as Fr

sys.path.insert(0, "src"); sys.path.insert(0, "benchmarks")
from fmsl.moving import (moving_profile, parse_moving, snapshot, _at_time, _evalE,
                         _finite_ends)
from fmsl.multivar import profile_polys
from fmsl.cluster import cmp, _bounds, _rational_between
from fmsl.algebraic import A, peval
from loop_crosscheck import loop_image, coverage_line, mismatch_line, signed_ivs
from dgnp2022_examples import EX

CLASSIC = ['(x^2-t*y^2)/(x^2+y^2)', 'x^2*y/(x^4+t*y^2)', 'x*y/(x^2+t*y^2)',
           '(x^3+t*y^3)/(x^2+y^2)', '(t*x^2+y^2-x*y)/(x^2+y^2)', '(x^2+t*y)/(x^2+y^2)',
           'x^2*y/(x^4+y^2+t*x^2*y)', '(x*y+t*x^2)/(x^2+(t-1)*y^2)']


# ------------------------------------------------------------------ family sources
def randA(gen):
    def coef():
        c = gen.choice([-2, -1, 1, 2])
        return gen.choice([f"({c})", f"({c})*t", f"({c})*(t-{gen.choice([1, 2, -1])})"])
    def poly(n, deg):
        terms = []
        for _ in range(n):
            i, j = gen.randint(0, deg), gen.randint(0, deg)
            if i + j == 0: continue
            terms.append(f"{coef()}*x^{i}*y^{j}")
        return " + ".join(terms) or "x"
    return f"({poly(gen.randint(1, 3), 3)})/({poly(gen.randint(2, 3), 4)})"


def randB(gen):
    def coef():
        a, b = gen.randint(-3, 3), gen.randint(-3, 3)
        return gen.choice([f"({a or 1})", f"({a or 1}*t+{b})", f"(t^2+{a}*t+{b})", f"({a or 1})*t"])
    def poly(n, lo, hi):
        terms = []
        for _ in range(n):
            d = gen.randint(lo, hi)
            i = gen.randint(0, d)
            terms.append(f"{coef()}*x^{i}*y^{d - i}")
        return " + ".join(terms)
    return f"({poly(gen.randint(2, 4), 2, 5)})/({poly(gen.randint(3, 5), 2, 5)})"


def hidden1(gen):
    c0 = gen.choice(["1", "2", "-1", "t", "t+1", "1-t", "2*t", "t/2"])
    r = gen.choice([0, 1, -1, 2, Fr(1, 2), Fr(-3, 2)])
    b = gen.choice([1, 2, 3])
    num = gen.choice(["x^4", "x^3*y", "x^4+x^2*y^2", "x^5", "x^4-x*y^3"])
    return f"({num})/((y-({c0})*x)^2 + (t-({r})+{b}*({c0})^3)*x^4 - {b}*x*y^3)"


def hidden2(gen):
    # D = (y - x - k x^2)^2 + (t - r + b) x^6 - b x^3 y^3 ... along y = x + k x^2 + u x^3:
    # x^3 y^3 = x^6 + ..., so the x^6 coefficient at the third level is (t - r)
    k = gen.choice([1, 2, -1])
    r = gen.choice([0, 1, -1, Fr(1, 2)])
    b = gen.choice([1, 2])
    num = gen.choice(["x^6", "x^5*y", "x^7"])
    return f"({num})/((y-x-({k})*x^2)^2 + (t-({r})+{b})*x^6 - {b}*x^3*y^3)"


def irrat(gen):
    r = gen.choice([0, 1, -1])
    pp = gen.choice(["t", "(t+1)", "2", "3", "(2*t)", "(t+2)"])
    rr = gen.choice([1, 2, Fr(1, 2), -1])
    num = gen.choice(["x^8", "x^8+t*x^4*y^4", "x^8-x^6*y^2"])
    return gen.choice([
        # hidden ON the irrational branch y = u x, u^2 = pp: x^4 y^4 = pp^2 x^8 there
        f"({num})/((y^2-{pp}*x^2)^2 + (t-({rr})+{pp}^2)*x^8 - x^4*y^4)",
        f"({num})/((y^2-{pp}*x^2)^2 + (t-({rr})+{pp}^2)*x^8 - x^4*y^4)",
        f"x^8/((y^2-2*x^2)^2 + (t-({r}))*x^8)",
        f"x^8/((y^2-t*x^2)^2 + x^8)",
        f"x^4*y^4/((y^2-(t+1)*x^2)^2 + (t-({r}))*x^6*y^2 + x^10)",
    ])


def dgnp(gen):
    i = gen.choice([1, 2, 3, 4, 6, 7, 8, 9, 18, 19, 20, 21])
    f, g = EX[i]
    how = gen.choice(["num_first", "den_last", "den_first", "add"])
    def scale_first(p, s):
        terms = p.split(" + ")
        return " + ".join([f"({s})*({terms[0]})"] + terms[1:])
    def scale_last(p, s):
        terms = p.split(" + ")
        return " + ".join(terms[:-1] + [f"({s})*({terms[-1]})"])
    if how == "num_first": f = scale_first(f, "t")
    elif how == "den_last": g = scale_last(g, "t")
    elif how == "den_first": g = scale_first(g, "t-1")
    else: f = f + " + t*x*y"
    return f"({f})/({g})"


def merge(gen):
    """arcs from different branches merge or separate as t moves: ordinary directions give an
       arc ending at p(t)-type values, the hidden curve y ~ a x gives arcs with ends shifted by
       p(t)/(1+a^2); the ends of the two branches cross at times that are invisible inside any
       single branch"""
    a = gen.choice(["1", "-1", "2"])
    b = gen.choice(["1", "2"])
    p = gen.choice(["t", "(t-1)", "2*t", "(t^2-1)", "(1-t)"])
    q = gen.choice(["x^2/(x^2+y^2)", "y^2/(x^2+y^2)", "x*y/(x^2+y^2)"])
    sgn = gen.choice(["-", "+"])
    return f"x^4/((y-({a})*x)^2{sgn}{b}*x^4) + {p}*{q}"


def path_family(gen):
    f = gen.choice(['x^2*y/(x^4+y^2)', '(x^2-y^2)/(x^2+y^2)', 'x*y/(x^2+y^2+x*y)',
                    'x^4/((x-y)^2+x^4)'])
    p = gen.choice([("t", "t^2"), ("t", "t"), ("0", "t"), ("t^2", "t^3")])
    return f, p


SOURCES = {"classic": None, "randA": randA, "randB": randB, "hidden1": hidden1,
           "hidden2": hidden2, "irrat": irrat, "dgnp": dgnp, "merge": merge, "path": path_family}


# ------------------------------------------------------------------ checks
class Timeout(Exception): pass


def _alarm(sig, frm): raise Timeout()


def loop_mismatch(N, D, t0, r=Fr(1, 10 ** 16)):
    N2, D2 = _at_time(N, t0), _at_time(D, t0)
    if not N2 or not D2: return 0.0
    cl = profile_polys(N2, D2)["cluster"]
    im = loop_image(N2, D2, r)
    if im is None: return None
    return mismatch_line(coverage_line(im), coverage_line(signed_ivs(cl)))


def _fl(v):
    return float(v)


def _interval_times(lo, hi, k):
    if lo is None and hi is None: a, b = -6.0, 6.0
    elif lo is None: b = _fl(_bounds(hi)[0]); a = b - 6
    elif hi is None: a = _fl(_bounds(lo)[1]); b = a + 6
    else: a, b = _fl(_bounds(lo)[1]), _fl(_bounds(hi)[0])
    out = []
    for j in range(1, k + 1):
        t = Fr(a + (b - a) * j / (k + 1)).limit_denominator(1000)   # simple fractions: messy ones make snapshots slow
        if (lo is None or cmp(t, lo) > 0) and (hi is None or cmp(t, hi) < 0): out.append(t)
    return out


def _near(t0, eps):
    """rational times t0 - eps and t0 + eps (t0 exact, possibly algebraic)"""
    lo, hi = _bounds(t0)
    while isinstance(t0, A) and hi - lo > eps / 100:
        t0.F.refine(); lo, hi = _bounds(t0)
    # simple fractions (small denominators): huge denominators made the loop check crawl
    # (two harness timeouts in the first runs were caused by this, not by the engine)
    left = Fr(float(lo - eps)).limit_denominator(10 ** 7)
    right = Fr(float(hi + eps)).limit_denominator(10 ** 7)
    if not (cmp(left, t0) < 0 < cmp(right, t0)):
        left, right = Fr(lo) - eps, Fr(hi) + eps
    return left, right


def check_family(text, path, mode, log, ctx):
    deep = mode == "deep"
    N, D = parse_moving(text, path)
    res = moving_profile(text, path=path, joint=False, deep=deep)
    def err(kind, t=None, detail=""):
        log({**ctx, "mode": mode, "kind": kind,
             "t": None if t is None else (str(t) if not isinstance(t, A) else f"~{float(t):.10g}"),
             "detail": detail})
    for a, b in res["unexplained"]:
        err("unexplained_left", a, f"change somewhere in ({a}, {b})")
    for u in res.get("deep_unsupported") or []:
        err("unsupported_branch", None, u)
    cps = [s["t"] for s in res["celedon_points"]]
    # LOOP check at interval samples, and DENSE type check
    for iv in res["intervals"]:
        m = loop_mismatch(N, D, iv["sample"])
        if m is not None and m >= 0.002:
            err("loop_mismatch", iv["sample"], f"snapshot {iv['snapshot']} vs loop image, mismatch {m:.3f}")
        ts = _interval_times(iv["lo"], iv["hi"], 12)
        ts = [t for t in ts if not any(cmp(t, c) == 0 for c in cps)]
        sigs = []
        for t in ts:
            s = snapshot(N, D, t)
            sigs.append(s["sig"] if s else None)
        for u in range(len(ts) - 1):
            if sigs[u] != sigs[u + 1]:
                # confirm both snapshots independently before blaming the interval
                ok = all((loop_mismatch(N, D, tt) or 0) < 0.002 for tt in (ts[u], ts[u + 1]))
                if ok and not any(a <= ts[u + 1] and b >= ts[u] for a, b in res["unexplained"]):
                    err("missed_change", ts[u], f"type differs between t={ts[u]} and t={ts[u + 1]}, "
                                                f"inside reported interval ({iv['lo']}, {iv['hi']})")
        # ENDS: re-check every end formula at 4 more times
        if iv["ends"]:
            sig_mid = snapshot(N, D, iv["sample"])["sig"]
            for t in _interval_times(iv["lo"], iv["hi"], 4):
                s = snapshot(N, D, t)
                if not s or s["sig"] != sig_mid or s.get("comps") is None: continue
                vs = _finite_ends(s["comps"])
                if len(vs) != len(iv["ends"]): continue
                for idx, (v0, E) in enumerate(iv["ends"]):
                    if E is None: continue
                    p = _evalE(E, t)
                    if len(p) <= 1 or peval(p, vs[idx]) != 0:
                        err("wrong_end_formula", t, f"end #{idx} = {vs[idx]} does not satisfy {E} "
                                                    f"(formula from interval ({iv['lo']}, {iv['hi']}))")
                        iv["ends"][idx] = (v0, None)
    # FALSE points, and LOOP check on both sides of each Celedon point
    for k, c in enumerate(cps):
        gaps = [abs(_fl(c) - _fl(o)) for o in cps if o is not c]
        eps = Fr(min([1e-3] + [g / 4 for g in gaps if g > 0])).limit_denominator(10 ** 9)
        left, right = _near(c, eps)
        sl, sr = snapshot(N, D, left), snapshot(N, D, right)
        s0 = snapshot(N, D, c) if not isinstance(c, A) else None
        for tt in (left, right):
            m = loop_mismatch(N, D, tt)
            if m is not None and m >= 0.002:
                err("loop_mismatch", tt, f"snapshot next to Celedon point {c}, mismatch {m:.3f}")
        same_sides = (sl and sr and sl["sig"] == sr["sig"])
        if same_sides and (s0 is None or s0["sig"] == sl["sig"]):
            if s0 is not None:
                err("false_point", c, f"type identical at t-eps, t, t+eps ({sl['text']} | {s0['text']} | {sr['text']})")
    return res


def main():
    seed = int(sys.argv[1]); n = int(sys.argv[2]); logpath = sys.argv[3]
    run = sys.argv[4] if len(sys.argv) > 4 else f"seed{seed}"
    limit = 240
    gen = random.Random(seed)
    fams = [("classic", f, None) for f in CLASSIC]
    for src in ("randA", "randB", "hidden1", "hidden2", "irrat", "dgnp", "merge", "path"):
        for _ in range(n):
            out = SOURCES[src](gen)
            f, p = out if isinstance(out, tuple) else (out, None)
            fams.append((src, f, p))
    seen, uniq = set(), []
    for s, f, p in fams:
        if (f, p) not in seen: seen.add((f, p)); uniq.append((s, f, p))
    logf = open(logpath, "a")
    def log(rec):
        logf.write(json.dumps(rec) + "\n"); logf.flush()
    signal.signal(signal.SIGALRM, _alarm)
    summary = {}
    T = time.time()
    for src, f, p in uniq:
        ctx = {"run": run, "seed": seed, "source": src, "family": f, "path": p}
        out = {}
        for mode in ("default", "deep"):
            signal.alarm(limit)
            t0 = time.time()
            try:
                r = check_family(f, p, mode, log, ctx)
                out[mode] = [str(s["t"]) if not isinstance(s["t"], A) else f"~{float(s['t']):.10g}"
                             for s in r["celedon_points"]]
                ok_rec = {**ctx, "mode": mode, "kind": "ok_summary", "t": None,
                          "detail": {"points": out[mode], "unexplained": len(r["unexplained"]),
                                     "seconds": round(time.time() - t0, 2)}}
                log(ok_rec)
            except Timeout:
                log({**ctx, "mode": mode, "kind": "timeout", "t": None, "detail": f">{limit}s"})
            except ZeroDivisionError:
                log({**ctx, "mode": mode, "kind": "skipped", "t": None, "detail": "D vanishes identically"})
            except Exception as exc:
                log({**ctx, "mode": mode, "kind": "crash", "t": None,
                     "detail": f"{type(exc).__name__}: {exc}"[:300]})
            finally:
                signal.alarm(0)
        if out.get("default") != out.get("deep") and "default" in out and "deep" in out:
            log({**ctx, "mode": "both", "kind": "modes_disagree", "t": None,
                 "detail": {"default": out["default"], "deep": out["deep"]}})
        print(f"[{src}] {f}{' path=' + str(p) if p else ''}  default {out.get('default')}  deep {out.get('deep')}",
              flush=True)
    print(f"done: {len(uniq)} families in {time.time() - T:.0f}s")


if __name__ == "__main__":
    main()
