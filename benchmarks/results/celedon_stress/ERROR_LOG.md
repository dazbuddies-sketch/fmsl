# Celedon points: the error log

*MJ Celedon. Every mistake found while stress-testing Celedon points, kept in one place so the
errors can be studied as a source of their own. That includes FMSL's mistakes and **our own
mistakes in testing**.*

Raw data: `error_log.jsonl` (1,361 records, one JSON line each: run, seed, source, family, mode,
kind, time t, detail). Verified catalogs: `catalog_*.md`. Tools: `benchmarks/celedon_stress.py`
(runs and logs) and `benchmarks/celedon_errors.py` (re-verifies every logged point and classifies it).

## How the tests were run

* **Sources.** Each seed has 57–60 families, from eight sources: the classic families, two
  different random generators (`randA`, `randB`), changes hidden one zoom level deep (`hidden1`)
  and two levels deep (`hidden2`), zoom curves that are not rational in t (`irrat`), the published
  examples of Dinh–Guo–Nguyễn–Phạm (2022) with a parameter inserted (`dgnp`), and moving paths
  (`path`).
* **Independent checks.** None of these use the engine's Newton-polygon zoom:
  * the exact image of f on a tiny square loop around the point;
  * 12 extra times inside every reported interval;
  * a before/at/after test at every reported point;
  * every end formula re-checked at 4 more times.
* **Runs.**
  * Seeds 101 and 202 before the fixes.
  * The same seeds after the fixes.
  * Seed 101 repeated.
  * Seed 303 with new families.
  * Seeds 202 and 303 again after the harness fixes.
* **Repeatability.** The two post-fix runs of seed 101 are identical: 175 of 175 records match
  (timings excluded).

## Results after the fixes (seeds 101, 202, 303; about 177 families)

| | default mode | deep mode |
|---|---|---|
| silent miss (real change, no answer, no warning) | 6 | **0** |
| miss with a "change near" warning | 37 | **0** |
| false Celedon point ("phantom") | 32 | **0** |
| disagreement with the loop image | 0 | 0 |
| wrong end formula | 0 | 0 |
| honest "cannot prove this branch" | – | irrational zoom curves only (`irrat`, some `dgnp`) |
| harness timeout (engine alone answers in 1.5–4.4 s) | 3 | 3 |

All default-mode errors came from the two hidden-zoom sources. On the classic, random, published
and moving-path families, the two modes agree everywhere.

## Error classes, causes and status

| # | error | where | cause | status |
|---|---|---|---|---|
| E1 | change hidden in the zoom missed, or found only approximately | default mode | candidates come from top-level Newton data only | **fixed in deep mode** (`fmsl.celedon`), proof in `CELEDON_DEEP.md` |
| E2 | false Celedon point next to a hidden change | default mode | the interval containing the hidden change differs from its neighbor, so the change is blamed on the nearest visible candidate | **fixed in deep mode** |
| E3 | wrong end formula ("v = 1" instead of v = 1/t) | default mode | a formula was accepted after a check at one sample time | **fixed in both modes**: confirmed at all three sample times |
| E4 | false Celedon point at an irrational time, e.g. t = ∛(1/2) | **both modes** | at an irrational time with an irrational zoom curve (y = √2·x), the 2-variable engine returns a partial snapshot (`exact = False`, e.g. {0} instead of [0, 0.28]); the moving layer ignored the flag | **fixed in both modes**: an unproven snapshot is no longer evidence of a change, and is labeled "(partial)". Deep-mode phantoms went 14 → 0 |
| E5 | snapshot cost spikes at times with large denominators | engine performance | rational roots are found by listing all divisors, which explodes for messy fractions | **open** (performance only, answers stay correct) |
| E6 | changes when arcs from different zoom branches merge | both modes | no exact condition for it yet | **open**: next step; not observed in these tests |
| E7 | zoom curves irrational in t (c₀ = √t) | deep mode | the method handles roots in Q(t) only | **reported honestly** as an unsupported branch, sampling as before |

## Our own mistakes during testing (fixed, and logged so they are not repeated)

| # | mistake | effect | fix |
|---|---|---|---|
| M1 | test generator missing parentheses: (2t)³ written as 2·t³ | families were not hidden exactly as designed (still valid tests) | parentheses added; reran |
| M2 | first version of the root finder said "all branches covered" when some roots were complex at the sample time (c² + t) | the completeness claim would have been wrong for such branches | count all deg_c roots; regression test added |
| M3 | checker judged irrational points by their two sides only | real changes that happen only *at* an irrational point (range collapsing to {∞} for one instant) were called phantoms | checker uses the exact algebraic value and a *proven* snapshot there |
| M4 | harness built test times with denominators up to 10⁶–10⁹ | slow snapshots and harness timeouts that looked like engine failures | simple fractions; the remaining timeouts re-judged directly (deep mode correct on both families) |

## What the errors teach (patterns)

1. **A displaced echo, not noise.** All 32 default-mode phantoms sit on a visible top-level
   candidate, and every family with a phantom (28 of 28) is also missing a real change. The real
   change happens where default mode cannot look, and its signal shows up at the nearest place it
   can look. The displacement ranged from 0.027 to 2.0 in t (median 1.0). So a phantom is a
   pointer: there is a hidden change near it. That alone could be turned into a warning in
   default mode.
2. **Partial information treated as complete.** E4 and M3 are the same shape of mistake in two
   places: an answer known to be incomplete was used as if it were complete. The fix in both
   places was to carry "proven / not proven" along with every answer.
3. **Spikes in cost, not in correctness.** E5 and M4: the cost of an exact answer spikes with the
   size of the numbers involved (the denominators), even when the mathematics is simple. Messy
   inputs make the method slower, not wrong.
4. **Single-point checks are fragile.** E3 (one sample time) and M2 (one sample time for roots)
   both came from checking at a single point. Every check that now matters is made at several
   points, or exactly.
