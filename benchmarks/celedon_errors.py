"""
Error catalog for the Celedon point stress test: re-verify every logged problem independently,
classify it, and look for patterns.

For each family where something was logged, every Celedon point reported by EITHER mode is
re-judged from scratch: the snapshot type just before and just after the point (each snapshot
confirmed by the square-loop image, which does not use the engine's zoom) decides whether the
point is a real change.  Then:
  real point missing from a mode  -> MISSED   (silent if that mode also left no warning)
  point reported but not real     -> PHANTOM  (and how far it sits from the nearest real point)

Run:  python benchmarks/celedon_errors.py LOG.jsonl [LOG2.jsonl ...] > catalog.md
"""
import json, sys
from collections import Counter, defaultdict
from fractions import Fraction as Fr

sys.path.insert(0, "src"); sys.path.insert(0, "benchmarks")
from fmsl.moving import parse_moving, snapshot, candidate_times
from celedon_stress import loop_mismatch


def as_time(s):
    return Fr(s) if not s.startswith("~") else Fr(float(s[1:])).limit_denominator(10 ** 12)


def confirmed_sig(N, D, t):
    s = snapshot(N, D, t)
    m = loop_mismatch(N, D, t)
    return (s["sig"] if s else None), (m is None or m < 0.002)


def judge(N, D, t, others, cache, exact=None):
    """is t a real change point?  returns (real, confirmed_by_loop).  `exact` is the exact
       (possibly algebraic) value of t when t is only known approximately: a change that happens
       only AT an irrational point is judged from the engine's snapshot there, but only if the
       engine marks that snapshot as proven (exact)."""
    if t in cache: return cache[t]
    gaps = [abs(float(t - o)) for o in others if o != t]
    eps = Fr(min([1e-4] + [g / 4 for g in gaps if g > 0])).limit_denominator(10 ** 12)
    (sl, cl), (sr, cr) = confirmed_sig(N, D, t - eps), confirmed_sig(N, D, t + eps)
    real = sl != sr
    conf = cl and cr
    if not real:
        s0 = snapshot(N, D, exact if exact is not None else t)
        if s0 is None or s0.get("exact", True):
            real = (s0["sig"] if s0 else None) != sl              # isolated change at the point itself
    cache[t] = (real, conf)
    return cache[t]


def main(paths):
    recs = [json.loads(l) for p in paths for l in open(p)]
    fams = defaultdict(list)
    for r in recs:
        fams[(r["source"], r["family"], json.dumps(r["path"]))].append(r)
    rows, kinds = [], Counter()
    by_source = defaultdict(Counter)
    for (src, f, pjs), rs in fams.items():
        path = json.loads(pjs)
        path = tuple(path) if path else None
        summ = {r["mode"]: r["detail"] for r in rs if r["kind"] == "ok_summary"}
        for r in rs:
            if r["kind"] in ("crash", "timeout", "loop_mismatch", "wrong_end_formula", "missed_change"):
                rows.append({"source": src, "family": f, "path": path, "mode": r["mode"],
                             "class": r["kind"].upper(), "t": r["t"], "detail": r["detail"]})
            if r["kind"] == "unsupported_branch":
                by_source[src]["deep: honest 'cannot prove' branch"] += 1
        if "default" not in summ or "deep" not in summ: continue
        N, D = parse_moving(f, path)
        pts = {m: [as_time(s) for s in summ[m]["points"]] for m in ("default", "deep")}
        allp = sorted(set(pts["default"]) | set(pts["deep"]))
        approx = set(s for m in ("default", "deep") for s in summ[m]["points"] if s.startswith("~"))
        exact_of = {}
        if approx:
            from fmsl.celedon import all_candidate_times
            ex = [c for c, _ in all_candidate_times(N, D)[0]]
            for s in approx:
                v = float(s[1:])
                m = [c for c in ex if abs(float(c) - v) < 1e-7 * max(1.0, abs(v))]
                if m: exact_of[as_time(s)] = m[0]
        cache = {}
        verdict = {t: judge(N, D, t, allp, cache, exact_of.get(t)) for t in allp}
        real = [t for t in allp if verdict[t][0]]
        top = [float(c) for c, _ in candidate_times(N, D)]
        for m in ("default", "deep"):
            mine = set(pts[m])
            for t in real:
                if t not in mine:
                    silent = summ[m]["unexplained"] == 0
                    rows.append({"source": src, "family": f, "path": path, "mode": m,
                                 "class": "MISSED_SILENT" if silent else "MISSED_FLAGGED",
                                 "t": str(t) if t.denominator < 10 ** 6 else f"~{float(t):.8g}",
                                 "detail": "real change not reported" + (" and no warning given" if silent else
                                                                        " (a 'change near' warning was given)"),
                                 "confirmed": verdict[t][1]})
            for t in mine:
                if not verdict[t][0]:
                    near = min(real, key=lambda u: abs(float(u - t))) if real else None
                    is_top = any(abs(float(t) - c) < 1e-9 for c in top)
                    rows.append({"source": src, "family": f, "path": path, "mode": m,
                                 "class": "PHANTOM",
                                 "t": str(t) if t.denominator < 10 ** 6 else f"~{float(t):.8g}",
                                 "detail": ("top-level candidate" if is_top else "deep candidate")
                                 + (f"; nearest real change at {float(near):.6g}, offset {float(t - near):+.6g}"
                                    if near is not None else "; no real change anywhere reported"),
                                 "confirmed": verdict[t][1]})
    for r in rows:
        kinds[(r["mode"], r["class"])] += 1
        by_source[r["source"]][f"{r['mode']}: {r['class']}"] += 1
    n_fams = len({(r["source"], r["family"], json.dumps(r["path"])) for r in recs})
    print(f"# Celedon point error catalog\n\n{n_fams} distinct families, logs: {', '.join(paths)}\n")
    print("## Counts by mode and class\n\n| mode | class | count |\n|---|---|---|")
    for (m, c), n in sorted(kinds.items()): print(f"| {m} | {c} | {n} |")
    print("\n## By source\n")
    for s in sorted(by_source):
        print(f"- **{s}**: " + ", ".join(f"{k} {v}" for k, v in sorted(by_source[s].items())))
    print("\n## Every error\n\n| source | mode | class | t | family | detail | loop-confirmed |\n|---|---|---|---|---|---|---|")
    for r in sorted(rows, key=lambda r: (r["class"], r["mode"], r["source"])):
        fam = r["family"] + (f" path={r['path']}" if r["path"] else "")
        print(f"| {r['source']} | {r['mode']} | {r['class']} | {r['t']} | `{fam}` | {r['detail']} | {r.get('confirmed', '')} |")
    json.dump(rows, open(paths[0] + ".catalog.json", "w"), indent=1, default=str)


if __name__ == "__main__":
    main(sys.argv[1:])
