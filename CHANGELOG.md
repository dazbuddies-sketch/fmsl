# Changelog

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
