# Changelog

## 0.32.0
- `fmsl.witness` zooms. At a rational point c0 where the initial denominator vanishes, it
  substitutes x_j = c0_j t^w_j and x_i = t^w_i (c0_i + y_i), then searches the new function with
  t > 0. This reaches curves with two or more terms, such as x = y sqrt(-z). It also searches
  weights up to 8 where the orders of P and Q tie. Every witness still has an explicit curve
  (`curve_point`).
- Part 1 (the answer is connected) is now used. Of the gaps between proven pieces, at most one
  lies outside the answer, and the result says so (`connected_note`).
- Check `benchmarks/witness_check.py` (`docs/theory/witness_check_output.txt`):
  - three variables, 50 functions: never outside the exact answer. Equal to it in 39, and
    equal up to gaps narrower than 1e-11 in 8 more (connectedness says at most one of those
    gaps is missing). 3 misses, down from 9;
  - four variables, 40 functions: 504,867 witnessed values re-checked along their curves by
    exact rational arithmetic, none wrong. 19 proven to be the whole circle, 33 proven to
    have no limit.

## 0.31.0
- Any number of variables: `fmsl.witness.inner_bound`. Every value it reports is a limit value
  of f at the origin, proven by an explicit curve x = (c1 t^w1, ..., cn t^wn) (weighted initial
  forms), and whole intervals are proven along segments of such curves (exact Sturm root
  isolation and the intermediate value theorem). It is an inner bound: values can be missing.
  When it covers the whole circle it is exact, and when it holds two values it proves the limit
  does not exist. Check `benchmarks/witness_check.py`:
  - three variables, 50 functions: the inner bound always lies inside the exact engine's answer,
    and equals it in 41. The 9 misses are whole-circle answers reached only along curves that
    need more than one term;
  - four variables, 40 functions with no exact answer to compare: 23,654 witnessed values, each
    re-checked by exact rational evaluation of f along its curve at t = 10^-20, 10^-40, 10^-80,
    none wrong. 18 proven to be the whole circle, 27 proven to have no limit.
- `benchmarks/part4_test.py`, Part 4 stress test (60 three-variable functions, half built to
  break the conjecture): 57 agree, 2 unfinished, 1 disagreement. The disagreement refutes the
  stronger statement that only the slices along the axis matter: there, the denominator's zero
  surface hugs the axis at order z^3. Details in `docs/theory/n-variables.md`.

## 0.30.0
- Three variables: clean indeterminacy lines (`fmsl.indeterminacy`). When numerator and
  denominator both vanish along a line through the origin with independent differentials there,
  every value is a limit value (proved in `docs/theory/n-variables.md`, Part 2); FMSL checks this
  exactly before zooming. Cross-check `benchmarks/indeterminacy_crosscheck.py`: of 300 random
  functions, the test fired on 70; the three-variable engine agreed on 69, never disagreed, and did
  not finish one, which the theorem now answers at once.
- `docs/theory/n-variables.md`: fluid states in n variables. Proved: the answer is one connected
  arc of the circle; clean indeterminacy points force the whole circle. Reduced: where the ends of
  a proper arc come from. Open, with evidence: non-clean indeterminacy branches as limits of slice
  fluid states (7 of 7 non-trivial cases agree with the independent engine; a first, simpler
  guess was refuted).

## 0.29.0
- Blow-up along curves for analytic quotients (`fmsl.analytic`, the curve method): when the
  denominator B is an exact polynomial whose zero curves the numerator A avoids, FMSL proves
  |A_K| + |B| >= c r^K and shows that A/B and A_K/B have the same limit values, +oo and -oo told
  apart (chordal distance -> 0; where either is large, A has the sign of A_K and B is the same).
  sin(x)/(x-y) -> [-oo, +oo], (1-cos(x))/(x^2-y^2) -> [-oo, 0] u [1/2, +oo],
  sin(x-y)/(x-y)^2 -> {-oo} u {+oo}.
- Cross-check `benchmarks/curve_crosscheck.py`: 76 functions, 73 proven, direct 320-digit
  evaluation never leaves a proven set, 70 of 70 damaged answers (a piece dropped, the wrong sign
  of infinity, a moved end) rejected.

## 0.28.0
- `limit_at_origin(..., time_limit=s)`: stop and report "not proven within s seconds".
- `fmsl.memory` remembers questions not decided within a time limit (timed out, or left
  undecided by the engine): asked again with the same or a smaller limit, FMSL says so at once; a
  larger limit tries again. On the learning benchmark, the 30 questions asked again in a new
  session took 0 s in total instead of 558 s.

## 0.27.0
- `fmsl.memory` (off by default): proven answers saved per FMSL version and returned at once in
  later sessions; a strategy learner (UCB1 bandit) that orders proof strategies by past reward
  (1/(1 + seconds) on success). Learning never decides an answer.
- Three variables with a distance: several starting directions for the squaring map, each with a
  5-second probe, then the best one without limit; lowest-degree direction first (found by the
  learner). On 30 new formulas (60 s cap): 15 proven instead of 8, in about the same total time
  (520 s against 511 s). Time slices respect a caller's own alarm, and use an exception engines
  cannot swallow.
- `benchmarks/learning_benchmark.py`.

## 0.26.0
- `fmsl.reward`: FMSL as a reward function and problem source for training AI models.
  `problems(n, families, seed)` gives prompts with proven answers (families: rational, analytic,
  distance, infinity, one_side); `score(problem, reply)` gives 1 for an exact answer and otherwise
  0.8 x the overlap of claimed and true sets (on the circle [-oo, +oo]), with what is missing and
  what is extra.
- `check_claim` and the MCP tool `fmsl_check_limit` check answers at infinity and from one side
  (`point: ["oo"]`, `side: "x:+,y:-"`).
- LimitBench-Beyond (`benchmarks/limitbench/limitbench_beyond.jsonl`): 200 problems in five
  families (sin/exp/log, common factors along lines, distances, at infinity, one side), each
  proven by FMSL and confirmed by direct 320-digit evaluation; `evaluate.py --set beyond`, and the
  model tools now take point and side.

## 0.25.0
- Limits at infinity in two and three variables (`fmsl.regions`, `limit_at(expr, "oo")`): the
  inversion q -> q/|q|^2 turns |q| -> oo into a limit at the origin; rational functions stay
  rational and a distance sqrt(Q) becomes sqrt(Q(u))/|u|^2.
- Limits from one side in several variables (`limit_at(expr, point, side={"x": "+", "y": "-"})`):
  half-planes, quadrants and octants, by x = a +- u^2.
- The MCP tool `fmsl_limit` takes `point: ["oo"]` and `side: "x:+,y:-"`.
- Cross-check `benchmarks/region_crosscheck.py`: direct 320-digit evaluation at points beyond
  1e36 (at infinity) or with the required signs; 65 limits, no sampled value outside a proven set,
  39 of 39 damaged answers rejected (43 sets filled by the samples; the rest are mostly
  whole-line answers, where sampling leaves gaps; 1 not proven).

## 0.24.0
- Common factors along lines, for sin, exp, log, ... in two variables: every rational line through
  the point on which the denominator vanishes is tested; the orders of numerator and denominator
  across it are bounded exactly by restricting the formula to the line in Q[t][w]/(w^J) (parts
  that cannot be restricted exactly keep a proven lower bound), and the common power is divided
  out of both Taylor series. sin(x+y)/(x+y) -> {1}, (1-cos(x+y))/(x+y)^2 -> {1/2}.
- The evaluator behind `fmsl.analytic` works over any ring (series, or restrictions to a line).
- Cross-check: 126 functions, now including 30 random ones with a common factor along a line;
  no sampled value outside a proven set, 117 sets filled, 32 of 32 damaged answers rejected. The
  checker's direct evaluation now uses 320 digits (80 digits lost tiny values such as
  log(1 + 1e-90) and briefly flagged two correct answers).

## 0.23.0
- Weighted distances: sqrt(Q) for any positive-definite quadratic form Q in the distances to the
  point (sqrt(x^2 + 2y^2), sqrt(x^2 + xy + y^2), sqrt(x^2 + 2y^2 + 3z^2), ...) that takes a rational
  square value. A Q-orthogonal rational basis starting at such a vector, and the map
  X1 = u^2 - D2 v^2 - D3 w^2, X2 = 2uv, X3 = 2uw, give sqrt(Q) = s (u^2 + D2 v^2 + D3 w^2). Forms
  with no rational square value (3x^2 + 3y^2) are refused with the reason.
- Distances written as powers are recognized precisely ((P)^(k/2)); a coefficient such as 1/2 no
  longer triggers that path.
- Fix: exact root finding when a polynomial's algebraic coefficients become rational after
  simplification (an AttributeError in `norm_poly`).
- Cross-check `benchmarks/form_crosscheck.py`: 48 formulas; for all 38 decided ones, three
  different starting vectors give identical answers and direct high-precision evaluation never leaves
  the proven set (37 sets filled by the samples, 1 whole-line answer with sampling gaps); 9 forms
  refused (no rational square value), 1 not proven within 90 s.

## 0.22.0
- sin, cos, tan, exp, log, sqrt, sinh, cosh, tanh, atan and asin in two and three variables
  (`fmsl.analytic`): exact Taylor polynomials A_Ka, B_Kb of numerator and denominator, and a
  proven Lojasiewicz bound |B_Kb| >= c r^L (r^(2L)/B_Kb^2 bounded, by the rational engine) with
  Ka >= 2L - b and Kb >= max(L, 2L - a), so that the remainders change no limit value; then the
  rational engine proves the answer. When the denominator vanishes along curves, the answer is
  "not proven", never a guess.
- `limit_at_origin`, `limit_at`, `check_claim` and the MCP tools accept these functions, also
  together with a distance; distances may be written as powers, e.g. (x^2+y^2+z^2)^(3/2).
- Square roots that are not distances, such as sqrt(1 + x), are left to `fmsl.analytic`.
- Cross-check `benchmarks/analytic_crosscheck.py`: direct 80-digit evaluation of the original
  formulas (no Taylor polynomials) at about 12,000 points within 1e-6 of the point; 84 functions,
  no sampled value outside a proven set, 76 sets filled by the samples (the other 8 are the whole
  line [-oo, +oo], where sampling leaves gaps), and 23 of 23 deliberately damaged answers rejected.

## 0.21.0
- Distances in three variables: `fmsl.distance` removes sqrt(x^2+y^2+z^2) (or a distance in two
  of three variables, or in two variables) exactly by the squaring map x = u^2 - v^2 - w^2,
  y = 2uv, z = 2uw, under which the distance is u^2 + v^2 + w^2. The map is onto and |q|^2 is the
  distance, so the fluid state is that of a rational function at 0, proven by the existing engines.
- `limit_at_origin`, `limit_at` (distance to any rational point), `check_claim` and the MCP tools
  use it in two and three variables.
- Cross-check `benchmarks/distance_crosscheck.py`: in two variables the squaring map and the polar
  blow-up of `fmsl.radial` (unrelated methods) must agree; in three variables, direction-only
  functions against numerical sphere ranges, and general functions across the three choices of
  leading variable, with sampled values inside every proven set.

## 0.20.0
- Distances: `fmsl.radial.radial_limit` gives the exact fluid state of functions of x, y and
  r = sqrt(x^2+y^2) at a point, by the rational polar blow-up (generic directions as an exact
  one-variable image, special directions by the two-variable engine, t = oo in a second chart).
- `fmsl.verify.limit_at_origin`, `check_claim` and the MCP tools accept sqrt(x^2+y^2) in two
  variables.
- Cross-check `benchmarks/radial_crosscheck.py`: 74 functions agree with exact circle images at
  radius 1e-12; the check compares isolated values as well as coverage, so moved points are caught.

## 0.19.0
- One variable: `log`, `sqrt` and rational powers `f^(p/q)`, as exact power series when the
  argument tends to 1 (powers: to c > 0 with c^(p/q) rational), via the recurrences
  L' f = f' and g' f = a f' g. Logarithmic singularities and branch points are reported, never
  guessed. Covers finance formulas such as log(1+r)/r, ((1+r)^(1/12) - 1)*12/r, annuity factors.
- Safe code generation: the radius estimate now also sees radius-1 series (log, sqrt), so the
  series zone stays inside the region of convergence.
- Several variables: fractional powers are rejected with a clear message.

## 0.18.0
- `fmsl.verify.limit_at(expr, point, side=None)`: limits at any rational point (by moving it to
  the origin), and one-sided limits in one variable (side "-" from below, "+" from above). The
  MCP tool `fmsl_limit` takes `point` and `side`.
- The non-properness documentation now credits the equivalence "bijection iff proper" to
  Hadamard's global inverse function theorem, and cites Jelonek (1993) and Pinchuk (1994).
- JOSS paper revised: the new features in one paragraph, verification of every result type,
  under 1,000 words.

## 0.17.0
- The non-properness (Jelonek) set of a polynomial map of the real plane (`fmsl.jelonek`,
  roadmap #6): Newton-Puiseux zooming along arcs to infinity; exact lines, half-lines and curves.
  Cross-check `benchmarks/jelonek_crosscheck.py` (elimination plus exact solving near each
  candidate, no arcs): 153 maps, 612 generic lines, 304 candidate points, no disagreements;
  446 witness arcs confirmed.
- MCP tool `fmsl_nonproperness`.

## 0.16.0
- `fmsl.sympy_bridge` (roadmap #5): `limit_set(expr, vars, point)` and
  `multivariate_limit(expr, vars, point)` for SymPy expressions in two or three variables, at any
  point with rational coordinates; exact endpoints returned as SymPy numbers (radicals when
  possible, otherwise CRootOf); `NoLimit` carries the set of limit values. `pip install
  "fmsl[sympy]"`. CI installs SymPy and the bridge tests may not be skipped there.

## 0.15.0
- Exact critical values when the fibre is nonlinear (roadmap #3): critical points over an
  algebraic first coordinate are now exact real algebraic numbers (roots of a polynomial with
  coefficients in that number field), and the value is computed in the compositum. Before, they
  were 40-digit approximations marked "approx."; on the 200 random three-variable functions
  this made 19 answers exact (same sets as before). The cube cross-check still matches 51/51.
  A flagged 40-digit fallback remains only beyond number-field degree 16.
- `fmsl.optimize`: exact global minimum / maximum of a two-variable function over the plane, and
  whether it is attained (and where) or only approached at infinity. Cross-check
  `benchmarks/optimize_crosscheck.py`: 67 functions, no disagreements.
- MCP tool `fmsl_extremum`.

## 0.14.0
- The rate (`fmsl.rate`, roadmap #4): for a two-variable rational function with a limit L at the
  origin, the exact best exponent theta with |f - L| <= C r^theta, and a curve along which it is
  attained (a Lojasiewicz-type exponent), by Newton-Puiseux zooming on the real roots of the
  denominator's initial forms. `rate_1d` in one variable. Cross-check
  `benchmarks/rate_crosscheck.py`: 139 functions, 36/36 closed forms, no sample point refutes a
  bound, every reported curve attains it.
- `fmsl.verify.rate_at_origin` and the MCP tool `fmsl_rate`.

## 0.13.0
- `fmsl.verify`: `limit_at_origin(expr)` (the proven set of limit values at the origin, in one,
  two or three variables, with witness curves in two variables) and `check_claim(expr, claim)`
  (is a claimed limit or set right; which values are missing, which are never approached).
- `fmsl-mcp`, an MCP server with the tools `fmsl_limit` and `fmsl_check_limit`, so AI assistants
  can compute and check limits instead of guessing (`pip install "fmsl[mcp]"`).
- LimitBench (491 problems) and LimitBench-Gaps (225 problems built around gaps in the set of
  limit values), with results for Claude Haiku 5.5 and Sonnet 5.5 alone, with FMSL available,
  and told to check their answers with FMSL (Sonnet 225/225, Haiku 224/225).

## 0.12.0
- The scaling map (`fmsl.scaling`), the rigorous form of DUS (Dimensional Unit Scaling) from MJ
  Celedon's FMS: for a two-variable 0/0 point, the regimes of scales y ~ c|x|^a with their
  dominant terms and behaviour, the balance scales (corners of the Newton polygons), a flag where
  a finer scale (zoom) decides, and which scale decides each end of the fluid state.

## 0.11.0
- One variable: `abs` and `sign`, analysed exactly on each side of 0 (near 0, f = c x^k (1 + o(1))
  fixes its sign on each side, so |f| = +-f and sign(f) = +-1 there).  An undecidable sign is
  reported, never guessed; code generation for abs/sign is refused for now.
- Experiment A rerun: FMSL 1120 correct, 0 wrong, 0 unresolved (was 1008 / 0 / 112); SymPy
  unchanged at 963 / 147 / 10.  A repeat run of 0.10.0 reproduced its tables exactly.

## 0.10.1
- The special times of `fmsl.moving` are named **Celedon points**, by MJ Celedon, after the
  Celedon Point of the Fluid Mathematics System (a transition between mathematical regimes):
  the exact times where the fluid state of a family changes character.  The result key
  `celedon_points` is added; `special` remains as an alias.

## 0.10.0
- Experiment A of MJ Celedon's FMS research brief (`benchmarks/experiment_a.py`, run in GitHub
  Actions against SymPy 1.13.3): 1,120 expressions in the brief's categories, ground truth by
  construction.  FMSL 1008 correct, 0 wrong, 112 unresolved (piecewise, unsupported); SymPy
  963 / 147 / 10; naive floats 557 / 279 / 60.  Several variables: FMSL exact in 224/224, the
  path method 67/224.
- One variable: essential factors are kept in products, so x*sin(1/x) -> 0 (squeeze theorem;
  it was reported as oscillating) and x^k exp(c/x^m) is decided on each side.
- Sums and powers of essential factors are reported as not decided (sin(1/x) - sin(1/x) was
  reported as oscillating).  The library functions Exp/Sin/Cos still raise Essential unless
  called with keep=True.

## 0.9.0
- Moving targets (`fmsl.moving`): for a family f(x, y, t), optionally with the 0/0 point moving
  along a path, the fluid state at each time, exact formulas E(t, v) = 0 for how its ends
  move, the exact special times where it changes character (rational or algebraic), and the
  fluid state of the moving target at each special time (three-variable engine).
- Cross-check `benchmarks/moving_crosscheck.py`: 75/75 snapshots agree with the square-loop
  image, 36/36 moving states with the cube image, no unreported change of type.
- Two variables: ranges of initial forms with algebraic coefficients are now exact (they were
  approximated to 40 digits and flagged); the approximation remains only as a flagged fallback.

## 0.8.0
- Three variables: exact zooming at irrational points (issue #2), e.g. along y = sqrt(2) x,
  z = sqrt(3) x for x^4/((y^2 - 2x^2)^2 + (z^2 - 3x^2)^2 + x^6) (fluid state [0, +oo]).
  On 47 functions built around irrational special points, complete answers rose from 34 to
  42 and UNPROVEN fell from 2 to 0; every changed or upgraded answer agrees with the cube
  cross-check.
- Exact real algebraic numbers: roots of polynomials whose coefficients are themselves
  algebraic (Sturm sequences over Q(alpha), norms), and arithmetic across two fields through a
  primitive element (common fields up to degree 16; beyond that a branch is reported
  unresolved, never guessed).
- Faster: one sign evaluation per interval refinement; cross-field equality decided exactly by
  the characteristic-polynomial test after a few refinements.  The 200 random functions of
  `benchmarks/space_random.py` take 35 s in total (slowest 26 s).

## 0.7.0
- Three variables: zoom along curves c^v0 d^v1 = u0 (issue #1). Denominators that vanish to
  high order along a curve such as yz = x^2 are now decided exactly; e.g.
  x^2yz/((yz - x^2)^2 + x^6) has fluid state [-1/4, +oo].
- Fix (also present in 0.6.0): a repeated root of the denominator's leading form was treated
  as a pole surface. 2x^3y^2z^2/((y^2z - x)^2 + x^4) was reported as [-oo, +oo]; its limit is 0.
  Found by the cube cross-check.
- Common monomial factors of numerator and denominator are cancelled first.
- Cube cross-check: now includes the fluid states at 0/0 points inside a face (it had missed
  the value 1/4 of -xyz/((yz - x)^2 + x^4 + z^4); the engine was right).
- New benchmark `benchmarks/space_random.py`: the engine decides 200 of 200 random functions.

## 0.6.0
- Highlights: all 21 examples of Dinh-Guo-Nguyen-Pham (2022) reproduced, including the curve
  case they left as future work; cross-checks agree 113/113 (2-D) and 51/51 (3-D).
- Three variables: `fmsl.space` computes exact fluid states of rational functions of x, y, z.
- Exact range of two-variable rational functions over a square or the plane (`fmsl.plane`).
- Signed ranges (+oo and -oo kept apart); all 21 examples of Dinh-Guo-Nguyen-Pham (2022)
  reproduced, including the curve case they left open.
- Independent cross-checks: loop (2-D) and cube (3-D) images; mutation tests.
- Faster exact arithmetic: cached integer Sturm sequences, continued-fraction rationals.

## 0.5.0
- Exact fluid-state arcs (cluster sets), resolution stopping rule, bifurcation stability.

## 0.4.0
- Exact real algebraic numbers, two-variable fluid states, bifurcations, sampling study.
