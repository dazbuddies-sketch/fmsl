"""LimitBench-Gaps: the harder companion set.

The first LimitBench runs (October 2026) showed that current Claude models find the values a limit
can take, but the few answers they got wrong all had the same flaw: they said f approaches every
value, when a whole range of values is in fact never approached (a gap). Proving a gap means
showing that no path at all reaches a value, which is exactly what an exact engine does.

So this set is made of functions whose set of limit values has a gap, plus about one quarter
decoys (the same kinds of functions but with no gap, so "there is a gap" is not a winning default).
Ground truth is checked the same way as LimitBench: FMSL's exact engine and the independent exact
image on small squares must agree, or the item is dropped.

Families:
  G1  hidden quadratic in t = y/x^k:   (x^2k + c y^2 + b x^k y)/(x^k y)  -> [-oo, b-2sqrt(c)] u [b+2sqrt(c), +oo]
  G2  ratio of shifted squares:         (y^2 + a x^2k)/(y^2 - b x^2k)    -> [-oo, -a/b] u [1, +oo]
  G3  random P(x^k, y)/Q(x^k, y), P and Q homogeneous of one degree (sometimes with extra terms of
      higher order), whose set has a gap
  D   decoys: built the same way, but whose set is one interval

Run:  python benchmarks/limitbench/generate_gaps.py  -> benchmarks/limitbench/limitbench_gaps.jsonl
"""
import json, os, random, signal, sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import generate
from generate import item_from_fmsl as _item              # same verification as LimitBench

# Gap functions are often undefined along an axis (x y in a denominator); the path method then
# simply uses the paths along which f is defined.
_paths = generate.path_values
generate.path_values = lambda N, D: [v for v in _paths(N, D) if v is not None] or [None]


class _Slow(Exception): pass


def _alarm(*_): raise _Slow()


def item_from_fmsl(f, *args, seconds=20):
    """The LimitBench item, or None if FMSL and the check take longer than `seconds` (the item is skipped)."""
    signal.signal(signal.SIGALRM, _alarm); signal.alarm(seconds)
    try: return _item(f, *args)
    except _Slow: return None
    finally: signal.alarm(0)


def has_gap(it):
    return it is not None and len(it["fluid_state_float"]) >= 2


def is_decoy(it):
    if it is None or len(it["fluid_state_float"]) != 1: return False
    lo, hi = it["fluid_state_float"][0]
    return lo != hi and not (lo == float("-inf") and hi == float("inf"))


def main():
    gen = random.Random(4099)
    items = []
    def add(it, tier, cat):
        if it is None or any(o["f"] == it["f"] for o in items): return False
        it["tier"], it["category"] = tier, cat
        items.append(it); return True
    # G1
    for k in (1, 2, 3, 4):
        for c in (1, 2, 3, 5):
            for b in (-3, 0, 2):
                f = f"(x^{2*k}+({c})*y^2+({b})*x^{k}*y)/(x^{k}*y)"
                it = item_from_fmsl(f, "", 0, f"t + c/t + b with t = y/x^{k}")
                if has_gap(it): add(it, 5, "gap: hidden quadratic")
    # G2
    for k in (1, 2, 3):
        for a in (1, 2, 3):
            for b in (1, 2, 4):
                f = f"(y^2+({a})*x^{2*k})/(y^2-({b})*x^{2*k})"
                it = item_from_fmsl(f, "", 0, f"(s+a)/(s-b) with s = (y/x^{k})^2")
                if has_gap(it): add(it, 5, "gap: shifted squares")
    # G3 and decoys: f = P(x^k, y)/Q(x^k, y) with P, Q homogeneous of the same degree d, so that the
    # values come from the one-variable function P(1, t)/Q(1, t), t = y/x^k; half of them get extra
    # terms of higher weighted order, which change nothing at the origin but hide the structure.
    def hom(d, k):
        T = [f"({gen.choice([-3, -2, -1, 1, 2, 3])})*x^{k * (d - j)}*y^{j}" for j in range(d + 1) if gen.random() < 0.7]
        return T or [f"y^{d}"]
    def extra(d, k):
        i, j = gen.randint(0, 4), gen.randint(0, 3)
        if k * j + i <= k * d: i += k * d - (k * j + i) + gen.randint(1, 2)    # weighted order > k d
        return f"({gen.choice([-2, -1, 1, 2])})*x^{i}*y^{j}"
    want_gap, want_decoy, tries = 110, 40, 0
    while tries < 6000:
        n_gap = sum(1 for it in items if it["category"] == "gap: random")
        n_dec = sum(1 for it in items if it["category"] == "decoy: no gap")
        if n_gap >= want_gap and n_dec >= want_decoy: break
        tries += 1
        if tries % 200 == 0: print(tries, "tries:", n_gap, "gap,", n_dec, "decoy", flush=True)
        k, d = gen.choice([1, 2, 2, 3, 3, 4]), gen.choice([2, 2, 3])
        P, Q = hom(d, k), hom(d, k)
        if gen.random() < 0.5: (P if gen.random() < 0.5 else Q).append(extra(d, k))
        f = f"({' + '.join(P)})/({' + '.join(Q)})"
        try:
            it = item_from_fmsl(f, "", 0, f"values of P(1,t)/Q(1,t), t = y/x^{k}")
        except (ZeroDivisionError, ValueError, ArithmeticError, RecursionError):
            continue
        if has_gap(it) and n_gap < want_gap: add(it, 6, "gap: random")
        elif is_decoy(it) and n_dec < want_decoy: add(it, 6, "decoy: no gap")
    for n, it in enumerate(items): it["id"] = f"LBG-{n + 1:04d}"
    with open(os.path.join(HERE, "limitbench_gaps.jsonl"), "w") as fh:
        for it in items: fh.write(json.dumps(it) + "\n")
    print(len(items), "items;", dict(Counter(it["category"] for it in items)), ";",
          sum(1 for it in items if it["path_test"]["trap"]), "path-test traps;", tries, "random tries")


if __name__ == "__main__":
    main()
