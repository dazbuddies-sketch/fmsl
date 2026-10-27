# FMSL — exact fluid states of 0/0 points

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23682015.svg)](https://doi.org/10.5281/zenodo.23682015)
[![tests](https://github.com/dazbuddies-sketch/fmsl/actions/workflows/tests.yml/badge.svg)](https://github.com/dazbuddies-sketch/fmsl/actions/workflows/tests.yml)

By **MJ Celedon**. FMSL computes exactly what a function does at a 0/0 point: the limit with a
proof when it exists, and otherwise the complete set of values it approaches (its *fluid
state*), for rational functions of two and three variables. It also generates proven-safe code
for one-variable formulas with 0/0 traps.

**Cite:** Celedon, MJ (2026). *FMSL: exact fluid states of 0/0 points.* Zenodo.
https://doi.org/10.5281/zenodo.23682015

## Proven-safe code for formulas with 0/0 traps

Formulas like `sin(t)/t`, `(1 - cos t)/t²`, `(1 - e^(-10r))/r` and the quantum-mechanics
"φ functions" `(e^(iz) - 1 - iz)/(iz)²` are everywhere in robotics, graphics, finance and
physics simulation. At zero they divide by zero, and **near** zero, ordinary
floating-point code silently returns garbage:

```python
>>> t = 1e-8
>>> (1 - math.cos(t)) / t**2     # true value: 0.5
0.0
```

FMSL is a small formula language and compiler that:

1. **Proves** what a formula does at the trouble point: removable value, pole (with its exact direction to infinity), or essential singularity.
2. **Proves** whether an expression is identically zero, or says honestly that it cannot.
3. **Emits safe Python or C code** that uses an exact series near the trap and the original formula elsewhere.
4. **Attaches a proof of the maximum error** in the danger zone near zero (a *Siletz bound*), then compiles the code and measures it as well.
5. **Provides robust matrix functions for quantum dynamics** (propagators, φ functions and control gradients) that handle degenerate energy levels without thresholds.

Real and complex formulas are supported. The core has no dependencies (pure Python 3.9+);
the matrix module uses numpy.

## Install

```
pip install fmsl            # core language and compiler
pip install fmsl[matrix]    # plus matrix functions (numpy)
```

or from source: `pip install .`

## Example

```
# robotics: SO(3) rotation coefficient
fluid B(t) = (1 - cos(t))/t^2

at t -> 0: B                       # removable: value 1/2 [proven]
check sin(t)^2 + cos(t)^2 - 1      # proven_zero
check exp(t^61) - 1                # proven_nonzero (order 61)
check exp(sin(t))*exp(-sin(t)) - 1 # unproven (zero through order 60; no certificate)
if B == 0.5                        # TYPE ERROR: '==' on fluid values is not allowed
emit c B                           # safe C code + proven error bound, compiled and validated
```

```
# quantum: time-evolution phase integral as the frequency w -> 0
context T = 3
fluid prop(w) = (exp(-T*i*w) - 1)/(-i*w)
at w -> 0: prop                    # removable: value 3 [proven]
check exp(i*t) - cos(t) - i*sin(t) # proven_zero  (Euler's formula)
```

Run the examples:

```
fmsl examples/realworld.fms
fmsl examples/quantum.fms
```

## Results

"Measured" is the worst relative error against exact values at about 190 points from
10⁻¹² to 1 on both sides of zero. "Proven" is the Siletz bound: a mathematical guarantee
for **every** input in the series zone near zero. Proofs are worst-case, so they are looser
than measurements, but they hold for all inputs, not just tested ones. 1e-16 is full
double precision.

**Robotics, graphics, finance**

| Formula | Field | Proven value at 0 | Naive code | FMSL measured | FMSL proven (zone) |
|---|---|---|---|---|---|
| (1 − cos t)/t² | Robotics (SO(3)) | 1/2 | 1.0 (100% wrong) | 1.2e-16 | ≤ 4.7e-15 (\|t\| < 1) |
| (t − sin t)/t³ | Robotics (SO(3)) | 1/6 | 1.0 (100% wrong) | 1.7e-16 | ≤ 4.4e-15 (\|t\| < 1) |
| (1 − t·sin t / (2(1 − cos t)))/t² | Robotics (SO(3) log map) | 1/12 | inf / NaN | 4.9e-16 | ≤ 1.5e-14 (\|t\| < 1) |
| (1 − e^(−10r))/r | Finance (annuity, r → 0) | 10 | 1.1e-05 | 2.0e-16 | ≤ 1.8e-13 (\|r\| < 0.215) |
| sin(2t/3)/sin t | Graphics (slerp) | 2/3 | fails only at t = 0 | 1.7e-16 | ≤ 3.2e-15 (\|t\| < 1e-4) |

**Quantum mechanics** (ħ = 1, evolution time T = 3)

| Formula | Use | Proven value at 0 | Naive code | FMSL measured | FMSL proven (zone) |
|---|---|---|---|---|---|
| (e^(−iTw) − 1)/(−iw) | Time-evolution phase integral | 3 | 4.7e-09 | 9.1e-17 | ≤ 2.8e-14 (\|w\| < 0.5) |
| φ₁ = (e^(iz) − 1)/(iz) | Schrödinger / Lindblad solvers | 1 | 5.0e-09 | 5.8e-17 | ≤ 1.7e-14 (\|z\| < 1) |
| φ₂ = (e^(iz) − 1 − iz)/(iz)² | Schrödinger / Lindblad solvers | 1/2 | 1.0 (100% wrong) | 1.1e-16 | ≤ 1.2e-14 (\|z\| < 1) |
| φ₃ = (e^(iz) − 1 − iz + z²/2)/(iz)³ | Schrödinger / Lindblad solvers | 1/6 | **8.1e+12** | 3.0e-16 | ≤ 1.1e-14 (\|z\| < 1) |
| (e^(−iTd) − 1 + iTd)/d² | 2nd-order perturbation theory | −9/2 | 1.0 (100% wrong) | 5.1e-17 | ≤ 9.6e-14 (\|d\| < 1) |
| sin²(Td/2)/(d/2)² | Fermi's golden rule | 9 | fails only at d = 0 | 0 | ≤ 3.9e-14 (\|d\| < 1) |
| (e^(−Tg) − 1 + Tg)/g² | Open-system decay kernel | 9/2 | **1.3e+07** | 1.6e-16 | ≤ 1.8e-13 (\|g\| < 1) |

Soundness check: with the series deliberately truncated (2 to 16 terms) on six test
formulas, the proven bound stayed above the actual error in all 30 cases.

## Matrix functions for quantum dynamics

`pip install fmsl[matrix]` (adds numpy). For a Hermitian Hamiltonian H:

```python
from fmsl import matrix as fm
U  = fm.propagator(H, t)            # exp(-iHt)
P2 = fm.phi(2, H, t)                # phi_2(-iHt), for exponential integrators
dU = fm.frechet(H, dH, t)           # derivative of exp(-iHt) w.r.t. H in direction dH
grads = fm.frechet_many(H, [dH1, dH2, ...], t)   # quantum-control gradients, one eigen-decomposition
```

These use the eigen-decomposition (the fast method) with FMSL's safe kernels, so exact and
near-degenerate energy levels and zero energies need no thresholds. The divided difference
(e^a − e^b)/(a − b) is evaluated as e^b·φ₁(a − b), which never divides by a small gap.

The scalar kernel φ_k(z) works for all complex z: an exact-coefficient series with a
**proven** error bound on |z| ≤ 1 (≤ 1.3e-13 for k = 1 to 4), and elsewhere a per-element
choice between the series and the standard recurrence by running error analysis. Measured
against exact values on 2,300 points (oscillating, decaying and growing regimes): at most
6.6e-16 relative error for φ₁ to φ₄. The naive formula is off by up to 3.4e+21 (φ₄).

**Accuracy** — relative error (Frobenius norm) against an independent 200-bit reference,
8×8 Hamiltonians built exactly with prescribed eigenvalues, t = 1:

| Case | Quantity | FMSL | Textbook spectral formula | scipy (block-matrix) |
|---|---|---|---|---|
| Exact degeneracy + zero energy | gradient | 8.4e-16 | 6.0e-03 (with equality guard) | 3.0e-16 |
| | φ₂(−iHt) | 4.3e-16 | **5.4e-01** (with zero guard) | 2.6e-16 |
| Near-degenerate (gaps 1e-9, 1e-12) | gradient | 1.2e-15 | 1.7e-06 (equality guard) | 2.9e-16 |
| | φ₃(−iHt) | 6.6e-16 | **5.3e-01** | 2.7e-16 |
| Weak coupling (all \|E\| ≤ 1e-6) | gradient | 1.3e-15 | 5.6e-11 (even with tolerance 1e-8) | 2.2e-16 |
| | φ₂(−iHt) | 8.6e-16 | **3.5e-01** | 9.6e-17 |
| Generic spectrum | gradient | 1.2e-15 | 1.2e-15 | 2.6e-16 |

The textbook spectral formula fails because floating-point eigenvalues of a degenerate
Hamiltonian are never exactly equal, and a zero energy is never exactly zero, so the usual
guards do not fire. scipy's block-matrix routines (`expm_frechet`, augmented `expm`) avoid
eigenvalues entirely and are also accurate, slightly more so than FMSL.

**Speed** — what FMSL adds over scipy is that the *fast* spectral method becomes as robust
as the slow one. Quantum-control gradients for m = 10 control directions sharing one
Hamiltonian (random Hermitian H):

| n | FMSL `frechet_many` | scipy `expm_frechet` × 10 | Speedup | Agreement |
|---|---|---|---|---|
| 32 | 1.7 ms | 6.8 ms | 4.1× | 2.1e-15 |
| 64 | 4.3 ms | 28.8 ms | 6.7× | 2.9e-15 |
| 128 | 24.6 ms | 179.5 ms | 7.3× | 5.2e-15 |

Reproduce with `python benchmarks/matrix_bench.py`. Timings are from one machine and will vary.

## Two variables: fluid states with direction

In two variables, a 0/0 point has no single state: the answer depends on the path of
approach. `fmsl.multivar` computes the full **state profile** of a rational function of x, y
at the origin, exactly:

```python
from fmsl.multivar import describe
print(describe("x^2*y/(x^4+y^2)"))
# verdict: NO LIMIT (proven)
#   because value depends on the path constant c along y = c*|x|^2: (c) / (1 + c^2)
# balance points (critical curvatures a): 2
```

Every straight line into the origin gives 0 here, yet along the parabolas y = c·x² the value is
c/(1 + c²). The engine finds such hidden paths from the **Newton polygon**: the state can only
change at *balance curvatures*, where competing terms of the numerator and denominator tie.
Along each path y = c·|x|^a the problem reduces to one variable and is solved exactly; where a
leading form vanishes, the engine zooms in (Newton–Puiseux) to find paths hidden deeper, e.g.
y = x + c·x² for x⁴/((x − y)² + x⁴). Irrational balance constants such as c = √2 are handled
with exact real algebraic numbers (`fmsl.algebraic`: Sturm-sequence zero and sign tests, never
floating point).

**When to stop zooming (resolution rule).** Zooming is only needed at a *repeated* root of a
leading form, or a root shared by the numerator and denominator forms. A simple root is
already settled by a theorem:

- a simple root of the denominator form that the numerator form does not share means there is
  a real curve where D = 0 and N ≠ 0 (implicit function theorem). So f is unbounded near it. If f
  is near 0 on both sides of that curve, it also takes **every** real value nearby
  (intermediate value theorem);
- a simple root of the numerator form, in a sector where f is unbounded, means f changes sign
  through 0 there, so again it takes every real value.

This handles curves whose Puiseux series never ends, such as the pole curve x = −y³/(1 − y).
On 200 random functions it leaves **0 unproven** (21 before exact algebraic numbers and this
rule), and all 200 profiles take about 1 second.

**Why a LIMIT verdict is a proof.** The engine examines curved paths, not every possible approach.
That is enough by the **curve selection lemma** (Bruhat–Cartan–Wallace, Milnor): for a function
built from polynomials, if f does not tend to L, then f already fails to tend to L along a
single Puiseux arc x = ±tᵠ, y = c·tᵖ + … Every such arc falls into one of the profile's cases:
a sector between balance curvatures, a balance curvature with c not a root, or a zoom at a
root. So when every case gives L, the limit is L.

### The fluid state is an exact arc

Cancel the common factors of N and D, and treat f as a map into the circle R ∪ {∞}. Near the
origin it is continuous, and a punctured disc is connected. So the set of all values f
approaches (its **cluster set**) is closed and **connected**. There are only three shapes:

| Shape | Meaning | Example |
|---|---|---|
| a single point | LIMIT = L, or DIVERGES (the point ∞) | x³/(x² + y²) → {0} |
| a closed arc | no limit; exactly these values are approached | x²y/(x⁴ + y²) → [−1/2, 1/2] |
| the whole circle | every value, including ∞, is approached | x³/(y − x²) |

The profile now outputs this set exactly (`profile(text)["cluster"]`, printed by `describe`).
Its endpoints can be irrational: for x²/(x² + 2y²) + y²/(3x² + y²) the set is [0.579796…, 1],
where the lower end is an exact algebraic number. An arc can also pass through ∞: for
xy/(x² − 3xy + y²) the values form the arc from −1/5 up through ∞ to −1. Equality of algebraic
numbers from different fields is decided exactly, using the characteristic polynomial of
multiplication by the number.

On the 200 random functions the shapes split 85 single points, 13 proper arcs and 102 whole
circles. So when a random rational function has no limit, it usually takes every value near
the origin.

**Two independent checks.**

1. *The theorem itself.* A computed cluster set that is not connected would contradict the
   theorem, so some family of paths would be missing. All 200 random cluster sets are connected
   and exact. To measure how much this check can catch, I deleted parts of correct profiles on
   purpose. It caught **72%** of deleted balance curvatures (29 of 40), but **0%** of deleted
   pole and zero curves (0 of 36). Losing those leaves a smaller set that is still connected.
   So it is a partial detector, not a completeness proof.
2. *A different algorithm.* `benchmarks/loop_crosscheck.py` computes the image of f on small
   squares |x|, |y| = r around the origin. On each side f is a one-variable ratio of
   polynomials, monotone between exactly isolated breakpoints, and values at critical points are
   evaluated as exact algebraic numbers. No Newton polygon is used. As r → 0 the image must
   converge to the cluster set. On 13 classic and 100 random functions, the
   loop image matches the engine's exact arc at every scale tested (r = 10⁻⁴, 10⁻⁸, 10⁻¹⁶, 10⁻³²):
   **113 of 113, no difference on the circle.** This includes a curve hidden three zoom levels
   deep, whose peak on the loop is narrower than r². When I injected errors on purpose, the loop
   check caught **42 of 42**: all 19 deleted pole and zero curves, and all 23 deleted balance
   curvatures. The connectedness test caught 14 of those 42.

**Validation of path values.** Every reported path value is re-checked by direct exact
evaluation of the original function along that path. The original validation (17 classic
cases, 200 random functions) found no disagreements across 20,360 path checks. The v0.4 engine
was re-validated on 200 new random functions, with 0 disagreements across 6,810 path checks; for
18 functions, exact re-evaluation at algebraic points exceeded the checker's 4-minute limit.
Development uncovered and fixed three real errors (missed denominator roots, missed
numerator–denominator balance points, and a thin hidden curve that random sampling never
hits); all three are regression tests.

`python benchmarks/fluid_state_figures.py` draws the pictures from the engine's output.

### Reproducing the published examples

`benchmarks/dgnp2022_examples.py` runs all 21 examples of Đinh, Guo, Nguyễn and Phạm (2022) through
FMSL. **All 21 agree** with the paper, including example 20, where the paper reports that
Maple 2021's `limit` returned a wrong range (Maple gave (3/4, 1); the true range is [3/4, +∞]).
The paper did not compute a range for example 5, because its denominator vanishes along curves
("future work"). FMSL gives **[−∞, −3] ∪ [3, +∞]**, which is disconnected in the extended line
but connected on the circle. The independent loop computation confirms it.

The signed range (+∞ and −∞ kept apart, as in the paper) is `profile(...)["cluster"]["signed_text"]`.
Checked independently, the signed ranges agree with the loop computation on **53 of 53**
functions: 13 of the paper's examples and 40 random ones. The 8 very-high-degree product
examples (10–17) were left out of the loop check for time; FMSL's results on them still match
the paper.

## Three variables

`fmsl.space` extends the fluid state to rational functions of x, y, z:

```python
from fmsl.space import describe3
print(describe3("(x*y+y*z+z*x)/(x^2+y^2+z^2)"))
# verdict: NO LIMIT (proven)
# fluid state (signed): [-1/2, 1]
```

An arc into the origin with x ≢ 0 has the form x = ±sᵠ, y = c·sᵖ + …, z = d·sʳ + …, so it
carries a weight (1, a, b). The (a, b) quarter-plane is cut by the lines where two terms tie
(the normal fan of the Newton polyhedra of N and D, plus N–D ties) into:

- **faces**, where the forms are monomials, so f → 0, ∞ (with signs) or a constant;
- **edges**, where the forms depend on one monomial u = cᵛ¹dᵛ², a one-variable problem solved
  exactly;
- **vertices**, where the forms are genuine polynomials in (c, d). The exact range of fN/fD over
  the whole plane is computed by `fmsl.plane`: critical values from resultants, values at
  infinity from charts analysed by the two-variable engine, then one exact test value per gap,
  decided by a cylindrical real-root test.

Simple roots are settled by the implicit function and intermediate value theorems. Repeated or
shared roots are zoomed: a change of variables, then the same analysis. Arcs inside the
coordinate planes go to the two-variable engine. The connectedness theorem holds in any number
of variables, because removing a curve of indeterminacy from a 3-ball leaves it connected. So
the result is again a point, an arc or the whole circle, and connectedness is checked.

Examples, each confirmed by hand or by a known closed form:

| f(x, y, z) | Fluid state | Why |
|---|---|---|
| xy/(x² + y² + z²) | [−1/2, 1/2] | 2xy ≤ x² + y² |
| x²y/(x⁴ + y² + z²) | [−1/2, 1/2] | the parabola trap, in 3-D |
| x⁴/((x − y)² + x⁴ + z²) | [0, 1] | curve hidden two levels deep |
| (xy + yz + zx)/(x² + y² + z²) | [−1/2, 1] | eigenvalues of the quadratic form |
| (xy − 2z² + xz)/(x² + 2y² + 3z²) | [−0.796873…, 0.394973…] | generalized eigenvalues |
| x²y²z²/(x⁶ + y⁶ + z⁶) | [0, 1/3] | decided at a vertex (c = d = 1) |
| xyz/(x⁴ + y⁴ + z²) | [−√2/4, √2/4] (to 40 digits) | max of cd/(1 + c⁴ + d²) |
| x²yz/((yz − x²)² + x⁶) | [−1/4, +∞] | zoom along the curve yz = x²: t/(t − 1)² ≥ −1/4 |
| 2x³y²z²/((y²z − x)² + x⁴) | {0} (limit 0) | a squared form touches 0 without crossing it |

On 200 random rational functions of three variables (`benchmarks/space_random.py`), the engine
decides **all 200** (24 limits, 162 no-limit, 14 divergent; median time 0.01 s). Every computed
fluid state is connected.

**Independent check.** `benchmarks/cube_crosscheck.py` computes the exact image of f on the
surface of small cubes |x|, |y|, |z| ≤ r. On each face f is a two-variable function over a
square, and its exact image comes from the box version of the range algorithm, plus the
two-variable fluid state at any point inside a face where the face function is 0/0. No Newton
polyhedra, weights or zooming are used. On 11 classic and 40 random functions, and on 54
functions built to vanish along curves (the curve zoom), the cube image matches the engine's
exact signed fluid state at r = 10⁻³ and 10⁻⁶ on **102 of 105**. The other 3 were too slow for
the cube (over 2 hours each); all three have the answer [−∞, +∞], confirmed independently: one
on the plane z = x, where the two-variable engine proves every real value is taken, and two by
explicit families of paths along which f tends to any chosen real number. The cube check found one real
engine error, fixed in 0.7.0: a squared denominator form had been treated as a pole surface.
When I deleted parts of correct answers on purpose (vertex ranges, pole and zero surfaces), the
cube check caught **11 of 11** changed answers. The connectedness test caught none of them. The low-level range algorithm itself was
tested against dense sampling: 60 of 60 random functions agree over a square, and 60 of 60 over
the whole plane.

Zooms at irrational points are exact: the analysis continues over the number field of the
point (for example Q(√2, √3) for a zoom at y = √2·x, z = √3·x), using exact real algebraic
numbers throughout. On 47 functions built around irrational special points, 42 are decided
completely; every answer that this changed agrees with the cube check below.

**Version 1 limits.** The following are reported as UNPROVEN, never guessed: zooming along a
curve c^p d^q = u₀ with neither |p| nor |q| equal to 1; number fields of degree above 16;
depth > 4. A
critical value whose critical point is not cut out linearly is computed to 40 digits and
marked "approx.".

## Iteration: bifurcations are balance points

Iterating a map x → f(x, r), new orbits are born at bifurcations. At a bifurcation point the
period polynomial Qₙ = (fⁿ(x) − x)/(f(x) − x) is a 0/0 point, and the new orbit grows along a
branch e = c·|u|^a (u = x − x₀, e = r − r₀) that `fmsl.bifurcation` finds exactly from the
Newton polygon. Verified against brute-force iteration (millions of steps):

| Bifurcation | Exact branch | Iteration ÷ prediction at e = 10⁻⁵ | Measured exponent |
|---|---|---|---|
| Transcritical (logistic, r = 1) | e = u | 0.99999 | 1.000 |
| Saddle-node (x + r − x²) | e = u² | 1.00000 | 0.500 |
| Pitchfork (x + rx − x³) | e = u² | 1.00000 | 0.500 |
| Period doubling (logistic, r = 3) | e = 9u² | 1.00000 | 0.500 |
| Second period doubling (logistic, r = 1 + √6) | e = c·u² | 0.99999 | 0.500 |

For the second period doubling the coefficient is exact:
**c = 20 + 10√2 + (35/3)√3 + (35/3)√6 ≈ 82.9268**, computed in Q(√2, √3) and confirmed by
iteration to five digits. These are leading-order laws: the agreement improves as e → 0.

**Stability, exactly.** Whether a new orbit attracts is decided by its multiplier
λ = (fⁿ)′(x) along the branch, which is a Puiseux series λ₀ + k·sᵐ + … with |λ₀| = 1 at the
bifurcation. `branch_stability` computes it exactly. Saddle-node: one branch stable, one unstable.
Pitchfork, transcritical and both period doublings: the new branch is stable. At r = 3 the
period-2 multiplier is exactly 1 − 36u² = 1 − 4e, which matches the known closed form
−r² + 2r + 4. At r = 1 + √6 the period-4 multiplier is exactly **1 − 8√6·e + …**, the value
flip-bifurcation theory predicts (4 × |d(multiplier of the 2-cycle)/dr| = 4 × 2√6). Iteration
at e = 10⁻⁴ measures 0.9980403, against the prediction 0.9980404.

`python benchmarks/bifurcation_demo.py` reproduces the table.

## Repetition vs exact analysis

What do you learn by asking the same question many times along random approach paths
(Monte Carlo, 1,000 paths per function), compared with the exact profile? On 300 random rational
functions, sampling matched the exact verdict 228 times (76%):

| Exact verdict | Sampling agrees | Sampling says instead |
|---|---|---|
| LIMIT (79) | 71 | NO LIMIT: 8 |
| NO LIMIT (179) | 124 | LIMIT: 43, DIVERGES: 12 |
| DIVERGES (42) | 33 | NO LIMIT: 9 |

Each kind of disagreement teaches something:

- **Thin curves are invisible to sampling.** When exact says NO LIMIT and sampling says LIMIT,
  the deciding behaviour sits on a curve that random paths almost never hit. Example:
  2x²y²/(3y² + x³) tends to 0 along every path sampling tries, but it blows up next to the
  curve 3y² = −x³. On the built-in traps, more repetitions do not help: for x⁴/((x − y)² + x⁴)
  and x(y − x³)²/((y − x³)² + x⁷), sampling finds nothing even with K = 10,000 paths.
- **Broad structure is found fast.** For the parabola trap x²y/(x⁴ + y²) and the √2 trap,
  detection reaches 100% by K = 30. For the parabola trap, the histogram of sampled values
  sketches the exact profile range [−1/2, 1/2].
- **Slow convergence looks like path dependence.** In all 8 cases where exact says LIMIT and
  sampling says NO LIMIT, the exact verdict is right. Along shallow curves (a ≈ 1/4), the
  deviation at sampling scale |x| ≈ 10⁻⁸ is still 8·10⁻², and it falls to 8·10⁻²⁵ only at
  |x| ≈ 10⁻¹⁰⁰.

So repetition is a good way to find candidate structure, but it cannot prove a limit or rule
out a hidden curve. The exact profile does both, and it also explains what the samples showed.
`python benchmarks/ensemble_vs_exact.py` reproduces the study.

## The language

| Statement | Meaning |
|---|---|
| `context NAME = EXPR` | Fixed exact parameter, e.g. `context T = 3` |
| `fluid NAME(VAR) = EXPR` | Define a value; the transition point is `VAR -> 0` |
| `at VAR -> 0: NAME` | Report the proven behavior at the point |
| `check EXPR` | Three-valued zero test: `proven_zero`, `proven_nonzero (order k)`, or `unproven` |
| `if NAME == ...` | Rejected by the type checker: equality on these values may be unprovable |
| `emit python NAME` / `emit c NAME` | Generate safe code with a proven error bound (C is compiled and validated if a C compiler is installed) |

Expressions: `+ - * /`, integer powers `^`, `sin cos exp`, exact decimal numbers, the
imaginary unit `i`, context parameters, and previously defined names, including
composition such as `B(t^2)`. Formulas containing `i` compile to complex code
(`cmath` in Python, C99 `double complex` in C).

Python API:

```python
from fmsl import x, Sin, Cos, Exp, I
print(((1 - Cos(x)) / x**2).describe())              # removable, value 1/2
print((Exp(I*x) - Cos(x) - I*Sin(x)).describe())     # C_e (certified: identically zero)
print((I / x).describe())                            # pole of order 1, heading 90 deg from the right ...
```

## How it works

Every value is an expression graph whose behavior at 0 is a **lazy Laurent series** with
exact coefficients in the Gaussian rationals Q(i), computed only as deep as a question needs.

**Zero certificate.** Certified values have the form E/Q, where Q is a polynomial and E is
built from `x^j · e^(q(x))` terms with polynomial exponents q (`sin` and `cos` use `±i·q`).
These building blocks are linearly independent, so a nonzero E must show a nonzero
coefficient within a computable window:

* **Linear exponents** (e^(λx), sin, cos, polynomials): E satisfies a constant-coefficient
  linear ODE of order N, so the window is N terms.
* **Polynomial exponents** (e^(x²), sin(x³), …): the window is found exactly by Gaussian
  elimination over Q(i).

If every coefficient in the window is zero, the value is **proven** identically zero.

**Siletz error bound.** For the emitted series branch, FMSL proves the maximum relative error:

1. A rigorous upper bound on |h(z)| over a complex circle |z| = r is computed through the
   expression graph (|a·b| ≤ |a|·|b|, |e^w| ≤ e^|w|, |sin w| ≤ sinh|w|, and dominant-term
   lower bounds for divisors, which also prove the divisor has no zeros inside).
2. Cauchy's estimate turns that into a bound on every series coefficient not yet computed,
   which bounds the truncation error.
3. A lower bound on |h| along the real interval is proven by subdivision, with exact
   evaluation at midpoints.
4. Floating-point rounding is bounded with Higham's Horner-scheme error analysis.

All bound arithmetic is exact rational arithmetic with outward rounding.

## Limits (by design, stated honestly)

* **Compositions** like `exp(sin x)` and **division by** `sin`/`cos`: coefficients are
  exact, but "is it identically zero?" is reported as `unproven`, never guessed. This
  boundary is fundamental: zero-testing general expressions is undecidable
  (Richardson's theorem).
* **Essential singularities** such as `exp(-1/t²)`, `sin(1/t)` and `exp(i/t)` are detected
  and explained, but no series code is emitted, because none exists.
* **Constants** must be exact (rational or Gaussian rational): `exp(1 + t)` would need e
  itself as a coefficient and is refused.
* **One variable**, transition point at 0; the variable is real (values may be complex).
* **Proof scope:** the Siletz bound covers the series zone near zero, where the danger is.
  The direct formula used farther out calls the platform math library, which has no
  formal accuracy guarantee, so that zone is measured, not proven.
* **Proven bounds are worst-case**, typically 10 to 1000× above the measured error.
  The compiler may add extra series terms to tighten a proof (for example 53 terms
  instead of 19 for the SO(3) log-map term).
* **Matrix module:** Hermitian Hamiltonians only (closed-system dynamics). The scalar φ
  kernel has a proof on |z| ≤ 1; matrix-level accuracy is measured, not proven, since it
  also depends on the eigen-decomposition. For non-Hermitian generators (e.g. Lindblad),
  use scipy's block-matrix routines.

## Prior work

The ingredients have a long history: lazy power series (McIlroy), zero-testing for
exponential polynomials via linear ODEs, limit computation for exp-log functions
(Gruntz's algorithm, used in SymPy and Maple), Cauchy estimates and rigorous ball/interval
arithmetic (e.g. the Arb library), and Higham's rounding-error analysis. The quantum φ
functions have their own numerical-analysis literature on stable evaluation; matrix
derivatives use the Daleckii–Krein formula, and scipy's `expm_frechet` implements the
Al-Mohy–Higham block method used as the baseline here. FMSL's
contribution is the combination: an always-defined treatment of 0/0 that never claims a
false zero, proven/unproven status enforced by the language, and automatic generation of
safe code that ships with a proof of its error.

For two-variable limits, the closest prior work is Đinh, Guo, Nguyễn and Phạm, *Limits of real
bivariate rational functions* (arXiv:2202.04889, 2022). It gives a Newton–Puiseux algorithm,
implemented in Maple, that computes the full range of limit values when the denominator has an
**isolated** zero; that range is a closed interval. For a denominator that vanishes along curves,
the authors note the range "may be not connected" in the extended line [−∞, +∞] and leave
computing it to future work. FMSL handles that case. Viewed on the circle R ∪ {∞}, where +∞ and −∞
are one point, the range is always connected. FMSL computes it exactly, as a point, an arc or the
whole circle, with algebraic endpoints, and checks it against an independent loop computation.
Related work computes limit existence for analytic quotients (arXiv:2102.01242) and Łojasiewicz
exponents (arXiv:2405.06302); the latter measure how fast a limit is approached. For three or more variables, the Extended Hensel Construction approach (core.ac.uk/works/9191392)
decides whether a limit exists, and computes it, when the denominator has an isolated zero. I did
not find a method that computes the full set of limit values in three variables, which
`fmsl.space` does (version 1, with the limits stated above). I have not done a full literature
search, so claims of novelty here are provisional.

## Credit

FMSL is by **MJ Celedon**, creator of the *Fluid Mathematics System* (FMS). It is the rigorous
continuation of FMS: each FMS goal that could be made rigorous appears here in a testable form.

| FMS idea | FMSL counterpart |
|---|---|
| Zero Intermedium (𝓩): division by zero as a transition state | Exact 0/0 states: removable value, pole, or essential singularity |
| Celedon Point (Cₑ): x − x as balance | `C_e`, returned only when *proven* identically zero |
| Infinity Flow (𝓘): infinity as directed flow | Poles report the exact direction to infinity from each side |
| Siletz Precision System (𝓢): structured precision | Siletz bounds: proven error enclosures for emitted code |
| Meta-Value Processor (𝓜): context-dependent values | `context` parameters |
| Uniqueness (𝓤): values keep their identity | Values carry their full construction (provenance graph) |

The original FMS formulation did not survive rigorous testing as a whole; the mathematics here
is the part of its goals that could be made precise, proven and independently cross-checked.

Development was assisted by an AI system (Anthropic's Claude), which wrote much of the code and
tests under the author's direction. Every result reported here is verified by tests and
independent cross-checks in this repository.

## License

MIT — free for any use. See [LICENSE](LICENSE).
