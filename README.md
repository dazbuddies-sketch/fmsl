# FMSL — exact fluid states of 0/0 points

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23682015.svg)](https://doi.org/10.5281/zenodo.23682015)
[![tests](https://github.com/dazbuddies-sketch/fmsl/actions/workflows/tests.yml/badge.svg)](https://github.com/dazbuddies-sketch/fmsl/actions/workflows/tests.yml)

By **MJ Celedon**. FMSL computes exactly what a function does at a 0/0 point: the limit with a
proof when it exists, and otherwise the complete set of values it approaches (its *fluid
state*), for rational functions of two and three variables. It also generates proven-safe code
for one-variable formulas with 0/0 traps.

**Try it in your browser, nothing to install:** [dazbuddies-sketch.github.io/fmsl](https://dazbuddies-sketch.github.io/fmsl/)
Type a formula and see every value it approaches, next to what the usual "try a few paths" test concludes.

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
pip install fmsl[mcp]       # plus the MCP server for AI assistants (fmsl-mcp)
pip install fmsl[sympy]     # plus the SymPy bridge (multivariable limits for SymPy expressions)
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
| xyz/(x⁴ + y⁴ + z²) | [−√2/4, √2/4] (exact) | max of cd/(1 + c⁴ + d²) |
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
depth > 4. Critical values whose critical points are not cut out linearly are exact too (since
0.15.0): the critical points are real algebraic numbers over the number field of their first
coordinate, and the value is computed in the compositum. Before, they were computed to 40
digits and marked "approx."; on the 200 random functions this made 19 answers exact, all with
the same sets as before. A 40-digit fallback remains only if the number field would exceed
degree 16, and is then flagged.

## Distances: limits with r = sqrt(x² + y²)

Distances are everywhere (potentials 1/r, norms, the Euclidean length of an error), and
r = sqrt(x² + y²) is not a rational function of x and y. FMSL handles it exactly with the polar
blow-up written in rational functions only, x = r(1 − t²)/(1 + t²), y = 2rt/(1 + t²), which
covers the punctured disc one-to-one. The fluid state is the union of a one-variable image
(generic directions) and the local fluid states at finitely many special directions, each proven
by the two-variable engine.

```python
from fmsl.radial import radial_limit
radial_limit("x*y/r")["set"]               # {0}
radial_limit("(r - x)/y")["set"]           # [-oo, +oo]   (tan of half the angle)
radial_limit("(r - x)/r")["set"]           # [0, 2]
radial_limit("x^2/(r^2 + r^3)")["set"]     # [0, 1]
radial_limit("x^2*y/r^3")["set"]           # [-2/(3 sqrt 3), 2/(3 sqrt 3)], exact algebraic ends

from fmsl.verify import limit_at_origin    # also takes sqrt(x^2+y^2) directly (and so does the MCP server)
limit_at_origin("a*b/sqrt(a^2+b^2)")["set"]   # {0}
```

### Three variables: sqrt(x² + y² + z²)

In `fmsl.verify` (and the MCP server) a distance is removed exactly by a squaring map, the real
form of z → z² and its three-dimensional cousin:

    x = u² − v² − w²,   y = 2uv,   z = 2uw,   and then   sqrt(x² + y² + z²) = u² + v² + w².

The map is onto and |q|² is the distance, so the values f approaches at the point are exactly the
values the rational function g(u, v, w) approaches at 0, which the three-variable engine proves.

```python
limit_at_origin("x*y*z/sqrt(x^2+y^2+z^2)^3")["set"]        # [-0.19245, 0.19245] = ±1/(3√3), exact
limit_at_origin("(x^2+y^2)/sqrt(x^2+y^2+z^2)")["set"]      # {0}
limit_at_origin("z/sqrt(x^2+y^2)")["set"]                  # [-oo, +oo]   (a cylinder: z left alone)
limit_at("(x-1)/sqrt((x-1)^2+(y-2)^2+z^2)", [1, 2, 0])     # other points: the distance to that point
```

Checked by `benchmarks/distance_crosscheck.py`: in two variables the squaring map and the polar
blow-up (unrelated methods) give identical sets on all 74 functions; in three variables, 36 of 36
direction-only functions match a numerical search over the sphere, and 13 of 13 general functions
give the same proven answer for all three choices of leading variable, with every sampled value
inside it (27 more were skipped because the engine could not finish a proof for some choice within
60 s; none disagreed).

### Weighted distances: sqrt(x² + 2y²) and other quadratic forms

The same idea removes the square root of any positive-definite quadratic form Q in the distances
to the point, as long as Q takes a rational square value s² at some rational vector v0. A rational
change of basis that starts with v0 and is orthogonal for Q gives Q = s²(X₁² + D₂X₂² + D₃X₃²), and

    X₁ = u² − D₂v² − D₃w²,   X₂ = 2uv,   X₃ = 2uw   give   sqrt(Q) = s·(u² + D₂v² + D₃w²).

```python
limit_at_origin("y/sqrt(x^2+2*y^2)")["set"]              # [-1/√2, 1/√2]
limit_at_origin("x/sqrt(x^2+x*y+y^2)")["set"]            # [-2/√3, 2/√3]
limit_at_origin("x*y*z/sqrt(x^2+2*y^2+3*z^2)^3")["set"]  # [-1/(9√2), 1/(9√2)]
```

When Q takes no rational square value (3x² + 3y², or 2x² + 3y², by Legendre's theorem), removing
the root would need an irrational constant, and FMSL refuses with that reason. Checked by
`benchmarks/form_crosscheck.py`: of 48 formulas (random forms in two and three variables), 38 were
decided; for each, three different starting vectors v0 give identical proven answers, and direct
320-digit evaluation of the original formula at about 12,000 points near the origin never leaves
the proven set. 9 forms were refused (no rational square value) and 1 was not proven within 90 s.

Square roots of quadratic forms are handled as above; `sqrt(1 + x)` and similar roots that are
not 0 at the point go to `fmsl.analytic` (below); anything else is refused with a clear message. Checked by
`benchmarks/radial_crosscheck.py`: 74 functions (14 hand-picked, 60 random) against the exact
image of the circle of radius 10⁻¹² — all agree, and deliberately damaged answers (an end moved
or a single value shifted) are all caught.

## sin, cos, exp, log in two and three variables

Most formulas in physics, engineering and finance are not fractions of polynomials. FMSL now takes
sin, cos, tan, exp, log, sqrt, sinh, cosh, tanh, atan and asin in two and three variables (with
arguments that vanish at the point, log(1 + u), and sqrt(c² + u)), and still proves its answer.

```python
from fmsl.verify import limit_at_origin
limit_at_origin("sin(x*y)/(x^2+y^2)")["set"]            # [-1/2, 1/2]
limit_at_origin("(1-cos(x*y))/(x^2+y^2)^2")["set"]      # [0, 1/8]
limit_at_origin("(x-sin(x))/sqrt(x^2+y^2)^3")["set"]    # [-1/6, 1/6]   (with a distance)
limit_at_origin("sin(x*y*z)/(x^2+y^2+z^2)^(3/2)")["set"]  # [-1/(3√3), 1/(3√3)]
```

How it is proven (`fmsl.analytic`). Write f = A/B with A, B analytic, of orders a and b at the
point, and take their exact Taylor polynomials A_Ka, B_Kb (all coefficients rational). If the
rational engine proves that r^(2L)/B_Kb² is bounded, that is |B_Kb| ≥ c·r^L, with Ka ≥ 2L − b and
Kb ≥ max(L, 2L − a), then

    |A/B − A_Ka/B_Kb| ≤ C·(r^(Ka+1+b) + r^(Kb+1+a)) / r^(2L) → 0,

so along every sequence the two quotients differ by something tending to 0, and they have exactly
the same limit values. The rational engine then proves the answer for A_Ka/B_Kb. This is the
Łojasiewicz route of Strzeboński (arXiv:2102.01242), here combined with FMSL's full sets of limit
values. When the denominator vanishes along curves through the point (sin(x)/(x − y)), no such L
exists and FMSL answers "not proven" instead of guessing.

Checked by `benchmarks/analytic_crosscheck.py`, which uses none of this: it evaluates each original
formula directly in 320-digit arithmetic at about 12,000 points within 10⁻⁶ of the point, on random
mixed scales and on curves y = c·x^k. Of 126 functions (36 hand-picked, 60 random, 30 random with a
common factor along a line), no sampled value lies outside the proven set; for 117 the samples fill
the set, and the other 9 are answers equal to the whole line [−∞, +∞], where random samples leave
gaps. All 32 deliberately damaged answers were rejected.

**Common factors along lines.** In sin(x + y)/(x + y) the denominator vanishes along a whole line,
but so does the numerator. In two variables FMSL looks for every rational line through the point
on which the denominator vanishes, and computes how often the numerator and the denominator
vanish across it: it restricts the formula to the line, exactly, in a ring of polynomials (t along
the line, w across it, truncated in w). Parts it cannot restrict exactly, such as sin(t) along the
line, keep only a proven lower bound on their order (F(s) = s·G(s) when F(0) = 0). The common power
is then divided out of both Taylor series exactly, and the method above applies:

```python
limit_at_origin("sin(x+y)/(x+y)")["set"]                    # {1}
limit_at_origin("(1-cos(x+y))/(x+y)^2")["set"]              # {1/2}
limit_at_origin("(x+y)*sin(x*y)/((x+y)*(x^2+y^2))")["set"]  # [-1/2, 1/2]
```

When the numerator does not vanish on the line (sin(x)/(x − y)), the values really do blow up
along it, and FMSL still answers "not proven" rather than guess.

## At infinity, and from one side

Two more questions reduce exactly to a limit at the origin, so every engine above answers them.

**At infinity.** The inversion q → q/|q|² maps the punctured plane (or space) onto itself with
|q/|q|²| = 1/|q|, so (x, y) → ∞ in any manner exactly when the new point tends to 0. Rational
functions stay rational, and a distance sqrt(Q) becomes sqrt(Q(u, v))/(u² + v²).

```python
from fmsl.verify import limit_at
limit_at("x*y/(x^2+y^2)", "oo")["set"]                # [-1/2, 1/2]
limit_at("(x^2+y)/(x^2+y^2)", "oo")["set"]            # [0, 1]
limit_at("sqrt(x^2+2*y^2)/(x+2*y^2)", "oo")["set"]    # [-oo, -1] u [0, +oo]
```

**From one side.** x → a from above is x = a + u² with u → 0 (from below, x = a − u²): onto the
half-line and preserving approach. Restricting several variables gives quadrants and octants.

```python
limit_at("x*y/(x^2+y^2)", [0, 0], side={"x": "+", "y": "+"})["set"]   # [0, 1/2]
limit_at("y/x", [0, 0], side={"x": "+", "y": "+"})["set"]             # [0, +oo]
```

The MCP tool takes `point: ["oo"]` and `side: "x:+,y:-"`. Functions such as sin at infinity, and
one-sided limits together with a distance, are refused for now. Checked by
`benchmarks/region_crosscheck.py`: direct 320-digit evaluation at points beyond 10³⁶ (at
infinity) or with the required signs (one side); of 65 limits, no sampled value lies outside the
proven set, and all 39 deliberately damaged answers are rejected. (A first run with points only
beyond 10⁶ flagged (x³y + xy³)/(x⁴ + y⁶): its limit at infinity is 0, but approached as slowly as
|q|^(−1/3), so the checker had to look much farther out.)

## Worked example: horizon or true singularity?

At the horizon of a black hole the usual formulas blow up. Is space itself singular there?
`examples/black_holes.py` asks FMSL (units G = c = 1, mass M = 1, so a non-spinning black hole
has its horizon at r = 2):

| formula | at | FMSL |
|---|---|---|
| metric component g_rr = 1/(1 − 2/r) | horizon r = 2 | +∞ from outside, −∞ from inside |
| metric component g_tt = −(1 − 2/r) | horizon r = 2 | 0 |
| curvature (Kretschmann scalar) 48/r⁶ | horizon r = 2 | **3/4: finite** |
| curvature 48/r⁶ | centre r = 0 | +∞ |
| spinning black hole (Kerr, spin a = 1): curvature in r and c = cos θ | the ring r = 0, θ = π/2 | **[−∞, +∞]** |

The metric blows up and flips sign at the horizon while the curvature, which does not depend
on coordinates, stays finite: the horizon is a coordinate singularity, the centre a true one.
At the ring of a spinning black hole the curvature has no limit at all: depending on the
direction of approach it tends to any value, which is a fluid state in the sense of this
package. FMSL decides the limits; what they mean physically is for the physics.

## Moving targets

What if the 0/0 point is not fixed? `fmsl.moving` takes a family f(x, y, t), optionally with
the target moving along a path (p₁(t), p₂(t)), and reports:

* **snapshots** -- the fluid state at each fixed time t, grouped into stretches of the same type;
* **how the ends move** -- an exact relation E(t, v) = 0 for each finite end v;
* **Celedon points** -- the exact times (rational or algebraic) where the snapshot changes
  character, and at each one the fluid state of the *moving* target: all limits as
  (x, y, t) → (0, 0, t₀) together, computed by the three-variable engine.

```
>>> from fmsl.moving import describe_moving
>>> print(describe_moving("(x^2-t*y^2)/(x^2+y^2)"))
snapshots (fluid state at a fixed time):
  -oo < t < -1:  at t = -3: [1, 3]
      ends move with t as:  v = 1;  v = -t
  -1 < t < +oo:  at t = -1/2: [1/2, 1];  at t = 2: [-2, 1]
      ends move with t as:  v = -t;  v = 1
Celedon points (the times where the snapshot changes character):
  t = -1:  snapshot {1};  moving target {1}
      why: a critical point meets a pole; numerator and denominator forms share a root
```

For x²y/(x⁴ + t·y²) the ends satisfy 4t·v² = 1 (that is, ±1/(2√t)) for t > 0, and every value
is taken for t ≤ 0. For (xy + t·x²)/(x² + (t − 1)y²) the Celedon point is the real root
t ≈ −0.754878 of t³ − t² + 1, where a gap in the range closes.

**What is proven.** Every reported fluid state is exact. Candidate Celedon points are the real
roots of exact polynomials in t: the coefficients of N and D, and the discriminants and
resultants of the initial forms (computed over Q[t]). Changes that arise only deeper in the
zoom are found by sampling three times in every interval and reported as "change near ...",
never hidden; this part is evidence, not proof. End formulas come from exact evaluation at
rational times and interpolation, with the degree bounded in advance and two extra times checked.

**Independent check** (`benchmarks/moving_crosscheck.py`, 8 classic and 30 random families):
snapshots agree with the square-loop image at sample times 75 of 75; moving states agree with
the cube image at rational Celedon points 36 of 36; dense sampling (12 times per interval)
found no unreported change of type. The closest prior work I found decides limits of
two-variable rational functions with Sturm's theorem (Zeng and Xiao, J. Symbolic Comput. 96, 2020); I did not
find a tool that computes exact fluid states for families with a parameter, but my search was
brief, so this is provisional.

## Polynomial maps: where properness fails (Jelonek set)

For a polynomial map F = (f, g) of the real plane, `fmsl.jelonek` computes its non-properness
set S_F exactly: the points F approaches along sequences that go to infinity. F is proper exactly
off S_F (Jelonek studied this set and showed how to compute it over the complex numbers
[Ann. Polon. Math. 58 (1993) 259–266]). A polynomial map of the plane whose Jacobian
determinant never vanishes is a bijection exactly when it is proper (Hadamard's global inverse
function theorem), that is, exactly when S_F is empty. This is the object behind the real
Jacobian conjecture: Pinchuk's counterexample [Math. Z. 217 (1994) 1–4] has a nonempty S_F.

```python
from fmsl.jelonek import jelonek_set
jelonek_set("x", "x*y")["text"]               # the line a = 0
jelonek_set("x*y", "x*y^2 - y")["text"]       # the line b = 0; the line a = 1
jelonek_set("x^2*y^2 + x*y", "y")["text"]     # the line b = 0, a in [-1/4, +oo)
jelonek_set("x*y", "x^2*y^2 + x")["text"]     # the curve (a, b) = (t, t^2), t real
jelonek_set("x + y^2", "y")["empty"]          # True: an automorphism
```
How: the same Newton-Puiseux zooming as at 0/0 points, now along arcs to infinity
x = ±1/s, y = c₁s^e₁ + c₂s^e₂ + … (and with x, y exchanged). f and g become Laurent polynomials in
s; F stays bounded only where every negative power cancels, which happens only at real roots of
initial forms, where FMSL zooms. Where both orders are ≥ 0, the limits are polynomials in the
free coefficient, giving exact pieces (A(t), B(t)) of S_F.

Checked by `benchmarks/jelonek_crosscheck.py`, which uses no arcs. On generic lines a = r and
b = r it finds every candidate point by elimination (where the number of solutions of F = (a, b)
can drop: leading coefficients of resultants), then decides each candidate by solving F = (a, b)
exactly near it and watching whether a real solution runs off to infinity. On 153 maps (6 with
known answers, 147 random), 612 lines and 304 candidate points (279 in S_F, 25 not), it agrees
with FMSL every time; and along all 446 witness arcs FMSL reports, |p| grows while F(p) tends
to the claimed point. The MCP server has the tool `fmsl_nonproperness`.

## Global minimum and maximum: is it attained?

`fmsl.optimize` finds the exact global minimum or maximum of a function of two variables over the
whole plane, and decides whether it is actually **attained** or only approached at infinity, a
question numerical optimizers cannot settle.

```python
from fmsl.optimize import global_min, global_max
global_min("x^2 + (x*y - 1)^2")   # infimum 0, NOT attained: approached only as (x, y) -> infinity
global_min("x^4 + y^4 - 2*x*y")   # minimum -1/2, attained at (x, y) = (-1/sqrt2, -1/sqrt2) and (1/sqrt2, 1/sqrt2)
global_min("x^2*y^4 + x^4*y^2 - 3*x^2*y^2 + 1")   # Motzkin: minimum 0 at (+-1, +-1)
global_max("x/(1 + x^2 + y^2)")   # maximum 1/2, attained at (1, 0)
```
How: the infimum m is the lower end of the exact range of f (critical values, values approached
at infinity, one exact test value per gap). If m is attained at a point, that point is critical,
so m is a critical value; conversely a critical value equal to m is attained there. So *attained
⇔ m is a critical value*, decided exactly. For polynomials, and rational functions whose
denominator never vanishes.

Checked by `benchmarks/optimize_crosscheck.py`, which uses neither critical values nor the range
algorithm: with the exact real-points test, f = a has no real solution for a just below m, f = b
has one for b just above m, and (for rational m) f = m has a solution exactly when FMSL says the
minimum is attained. On 67 functions (7 classics with known answers and 60 random
polynomials of degree 4, of which 23 do not attain their infimum) there are no disagreements.

## The rate: how fast a limit is approached

Knowing that f → L is weaker than knowing how fast. `fmsl.rate` computes the best exponent θ in

    |f(x, y) − L| ≤ C·r^θ   near the origin   (r = distance to the origin),

exactly, together with a curve along which the difference is of exact order r^θ, so no larger
exponent works. This is a Łojasiewicz-type exponent; it decides, for example, how fast an
iteration or a controller settles near the point.

```python
from fmsl.rate import rate
rate("x^2*y/(x^2+y^2)")["exponent"]        # 1
rate("x^2*y^2/(x^2+y^4)")["exponent"]      # 2, attained along x = c y^2
rate("x^7/((y-x^2)^2+x^6)")["exponent"]    # 1, attained only along y = x^2 + c x^3
```
The last example is the reason for exactness: away from the curve y = x² the difference is of
order r³, and no grid of sample points finds the thin region where it is of order r.

How: along curves y = c|x|^b the order of f − L is a piecewise linear function of b for all but
finitely many c, minimal at a corner of the Newton polygons; the exceptional c are real roots of
the denominator's initial form, where FMSL zooms (y = |x|^b (c + y₁)) and repeats. Roots of the
numerator's initial form only make f − L smaller, so they cannot lower the minimum.

Checked by `benchmarks/rate_crosscheck.py` on 139 functions (the 103 two-variable LimitBench
functions with a limit, and 36 that need zooming, with a closed-form answer): all 36 match the
closed form; exact rational sample points very close to the origin never contradict the bound;
and on every reported curve the measured order agrees with θ (to 0.02). One variable:
`rate_1d`; the MCP server has the tool `fmsl_rate`. Three variables are not covered yet.

## Experiment A: a head-to-head test from the FMS research brief

MJ Celedon's original FMS brief asks for "a benchmark of at least 1,000 mathematical expressions"
across ten kinds of local behaviour, compared with established methods "using independently
established ground truth", counting correct and false classifications and unresolved cases.
`benchmarks/experiment_a.py` does this with **1,120 expressions** (112 per category; the
parameter-dependent and multivariable categories are two-and three-variable problems). Ground
truth comes from how each expression is built (known series, AM-GM bounds, eigenvalues), never
from a method under test. It runs in GitHub Actions against SymPy 1.13.3
(`.github/workflows/experiment_a.yml`); the full tables are in `benchmarks/results/`.

| method | correct | wrong | unresolved |
|---|---|---|---|
| FMSL | 1120 | 0 | 0 |
| SymPy (one variable: `limit` from each side; several variables: limits along 8 standard paths) | 963 | 147 | 10 |
| naive floating point (probing at 10⁻⁴ and 10⁻⁷) | 557 | 279 | 60 (of 896) |

What this shows, and what it does not:

* **One variable: a tie.** FMSL and SymPy are both correct on all 896 one-variable cases:
  removable, poles, one-sided, oscillatory, cancellation, ill-conditioned, essential, and
  piecewise (|x|, sign). FMSL handles |f| and sign(f) exactly by analysing each side of 0
  separately (near 0, f ≈ c·xᵏ fixes the sign on each side); this was added in 0.11.0, before
  which the 112 piecewise cases were unsupported. The run is reproducible: two runs of the same
  version gave identical tables.
* **Several variables and parameters: the difference.** SymPy has no multivariable limit
  routine, so the baseline is the common practice of limits along paths. It found the exact set
  of limit values in 67 of 224 cases, and got the yes/no question "does the limit exist?" wrong
  or unresolved in 19, e.g. x³y/(x⁶ + y²): all eight paths give 0, but along y = x³ the value is
  1/2 (the fluid state is [−1/2, 1/2]). FMSL found the exact set in 224 of 224.
* **Floating point** fails on cancellation exactly as expected; the brief rightly warns that
  beating it proves nothing on its own, which is why SymPy is the main baseline.
* Running the experiment exposed one FMSL error before release: x·sin(1/x) was reported as
  oscillating (the limit is 0, by the squeeze theorem). It is fixed in 0.10.0, together with a
  related problem (sums of essential factors now say "not decided" instead of borrowing a verdict).

## LimitBench: AI models against proven answers

[`benchmarks/limitbench/`](benchmarks/limitbench/) asks for the complete set of values a function
of two or three variables approaches at the origin, with every answer proven by FMSL and confirmed
by the independent exact check. Models are run through the Anthropic API in GitHub Actions
(`.github/workflows/limitbench.yml`); every answer is saved and the IDs of all wrong answers are
listed in [`results/REGRADED.md`](benchmarks/limitbench/results/REGRADED.md). Exact answers,
October 2026:

| method | LimitBench (491) | LimitBench-Gaps (225) |
|---|---|---|
| FMSL | 491 | 225 |
| Claude Sonnet 5.5 | 489 | 218 |
| Claude Haiku 5.5 | 489 | 216 |
| limits along 8 standard paths | 104 of 450 two-variable items | 0 |

What this shows, and what it does not:

* **Current models are good at this.** On the main set both models were right about whether the
  limit exists every time, including all 132 problems where every standard path gives the same
  wrong impression. The claim "AI gets fooled where FMSL does not" is not supported.
* **Their mistakes have a shape.** Proving which values are *never* approached is the hard half.
  LimitBench-Gaps is built around it: 185 functions whose set of limit values has a gap, and 40
  decoys without one. All 7 of Sonnet's mistakes there were the same: it said every value is
  approached when a whole range is not. For example, −2x²/(x² + 2xy − 2y² + x³y) never comes
  near any value strictly between −4/3 and 0. Haiku made 2 mistakes of that kind and 7 others
  (wrong endpoints, a stray point, a gap in the wrong place, one reply with no answer).
* **The practical point.** A model is right about 96–97% of the time on these problems, and its
  answer has to be checked to know which time it is. FMSL is right every time, proves it, and
  is the tool that does the checking.
* The grader was fixed twice after the runs to read answers written as formulas (such as
  (1−√10)/2 and 2√3); all saved answers were regraded with `regrade.py`. Model outputs vary
  between runs; these are single runs.

## For SymPy users

SymPy's `limit` handles one variable. `fmsl.sympy_bridge` adds two and three variables, with
SymPy expressions in and SymPy objects out, and the complete set of limit values when the limit
does not exist:

```python
import sympy as sp
from fmsl.sympy_bridge import limit_set, multivariate_limit   # pip install "fmsl[sympy]"
x, y = sp.symbols("x y", real=True)
limit_set(x**3*y/(x**6 + y**2), (x, y))           # Interval(-1/2, 1/2)
limit_set((x**2 + y**2)/(x*y), (x, y))           # Union({-oo, oo}, Interval(-oo, -2), Interval(2, oo))
multivariate_limit((x - 1)*(y - 2)**2/((x - 1)**2 + (y - 2)**2), (x, y), (1, 2))   # 0
```
Endpoints are exact, as radicals when possible: example 4 of Dinh, Guo, Nguyen and Pham (2022)
comes back as `Interval(-sqrt(2)/4, sqrt(2)/4)`. Any point with rational coordinates works.
When the limit does not exist, `multivariate_limit` raises `NoLimit`, whose `.values` is the
set. Tested in CI with SymPy 1.14.

## FMSL for AI assistants (MCP server)

`fmsl-mcp` gives an AI assistant five tools over the Model Context Protocol:

* `fmsl_limit(expression)`: the complete set of values the expression approaches at the origin,
  proven, with the curves along which the end values are approached (two variables).
* `fmsl_rate(expression)`: how fast the limit is approached (see *The rate* above).
* `fmsl_extremum(expression, kind)`: the exact global minimum or maximum over the plane, and
  whether it is attained.
* `fmsl_nonproperness(f, g)`: the Jelonek set of a polynomial map of the plane.
* `fmsl_check_limit(expression, claim)`: checks the assistant's own answer, and says exactly what
  is missing and what is never approached, e.g. *Incorrect: the proven set is [-oo, -4/3] u
  [0, +oo]; these claimed values are never approached: (-4/3, 0).*

```
pip install "fmsl[mcp]"
```
Claude Desktop or Claude Code configuration:
```json
{"mcpServers": {"fmsl": {"command": "fmsl-mcp"}}}
```
The same checks are available from Python as `fmsl.verify.limit_at_origin` and
`fmsl.verify.check_claim`, for graders and test suites. Other points work by moving them to the
origin (x = 2 + h, h -> 0), which `fmsl.verify.limit_at` does for you, with one-sided limits in
one variable: `limit_at("1/(x-2)", 2, side="-")` gives -oo (from below), `side="+"` gives +oo.
`fmsl_limit` takes the same `point` and `side`.

**Does it help? Only when the model is told to check.** LimitBench-Gaps (225 problems), exact
answers:

| model | alone | FMSL available | told to check its answer with FMSL |
|---|---|---|---|
| Claude Sonnet 5.5 | 218 | 220 (FMSL called 75 times) | **225** (every problem checked) |
| Claude Haiku 5.5 | 216 | 217 (FMSL called 4 times) | **224** (224 of 225 problems checked) |

Given the tool, the models mostly answered without it: they are confident about limits. With
one instruction in the system prompt (*"Before you give your final answer, check it with
fmsl_check_limit. If FMSL says it is wrong, correct it."*) Sonnet made no mistakes, and Haiku's
one miss was a reply that ended without an answer and never called the tool. So the recommended
setup is to instruct the assistant to verify. Single runs; details and per-problem tool calls in
`benchmarks/limitbench/results/`.

## The scaling map (DUS)

MJ Celedon's FMS proposed *Dimensional Unit Scaling* (DUS): transformations across numerical
scales. Its rigorous form is the step at the heart of FMSL's engine: probe the 0/0 point at every
scale y ~ c·|x|^a (y has "dimension" a relative to x). At each scale the smallest-power terms
dominate, so f behaves like |x|^δ times a constant. `fmsl.scaling` shows this as a map:

```
>>> from fmsl.scaling import describe_scaling
>>> print(describe_scaling("(x^2-y^2)/(x^2+y^2)"))
scaling map (y ~ c|x|^a; size of f ~ |x|^delta):
  0 < a < 1        N ~ -y^2,  D ~ y^2:  f -> -1
  a = 1            BALANCE  N ~ -y^2 + x^2,  D ~ y^2 + x^2:  balance: f of order 1, depends on c
  a > 1            N ~ x^2,  D ~ x^2:  f -> 1
fluid state: [-1, 1]
```

* **Regimes** are ranges of a where one term of N and one of D dominate; f → 0, → ∞, or a value.
* **Balance scales** are the finitely many a where terms tie (corners of the Newton polygons).
  Only there can the leading terms interact, so that is where the interesting values come from.
  When the denominator's leading form vanishes at some c, the map says so: a finer scale (a
  zoom) decides there, as for x⁴/((x − y)² + x⁴), whose end value 1 lives on y = x exactly.
* **Which scale decides each end** of the fluid state is listed (the paths from the exact analysis).

The regimes are computed exactly from the exponents; this is weighted (quasi-homogeneous)
scaling, the standard idea behind Newton polygons and tropical geometry, shown to the user.

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

Expressions: `+ - * /`, integer powers `^` and rational powers `^(p/q)`, `sin cos exp log sqrt`,
`abs sign` (analysed on each side of 0), exact decimal numbers, the
imaginary unit `i`, context parameters, and previously defined names, including
composition such as `B(t^2)`. Formulas containing `i` compile to complex code
(`cmath` in Python, C99 `double complex` in C). `log f`, `sqrt f` and `f^(p/q)` are exact power
series when f tends to 1 (for powers, to any c > 0 with c^(p/q) rational), which covers the usual
finance formulas: `log(1+r)/r -> 1`, `((1+r)^(1/12) - 1)*12/r -> 1`, `(1 - (1+r)^(-12))/r -> 12`.
`log f` with f -> 0 and fractional powers of a vanishing quantity are reported, not guessed: they
are not power series.

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
* **Several variables with sin, exp, log, ...:** proven only when the denominator's zero at the
  point is isolated, or becomes isolated after cancelling common factors along rational lines
  (two variables). Otherwise, as for sin(x)/(x − y), (sin x − sin y)/(x − y) or a curve such as
  y = x², FMSL says "not proven". Arguments of sin, cos, exp, ... must be 0 at the point, of log
  equal to 1, and of sqrt a positive rational square; anything else is refused, since it
  would need constants such as sin 1 or e.
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
| Celedon Point (Cₑ): a transition between mathematical regimes | **Celedon points** of a family (`fmsl.moving`): the exact times where the fluid state changes character |
| Celedon Point (Cₑ): x − x as balance | `C_e`, returned only when *proven* identically zero |
| Infinity Flow (𝓘): infinity as directed flow | Poles report the exact direction to infinity from each side |
| Siletz Precision System (𝓢): structured precision | Siletz bounds: proven error enclosures for emitted code |
| Dimensional Unit Scaling (DUS): transformations across scales | The scaling map (`fmsl.scaling`): which scale y ~ c\|x\|^a decides each part of a fluid state |
| Meta-Value Processor (𝓜): context-dependent values | `context` parameters |
| Uniqueness (𝓤): values keep their identity | Values carry their full construction (provenance graph) |

The original FMS formulation did not survive rigorous testing as a whole; the mathematics here
is the part of its goals that could be made precise, proven and independently cross-checked.

Development was assisted by an AI system (Anthropic's Claude), which wrote much of the code and
tests under the author's direction. Every result reported here is verified by tests and
independent cross-checks in this repository.

## License

MIT — free for any use. See [LICENSE](LICENSE).
