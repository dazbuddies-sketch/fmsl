# Changelog

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
