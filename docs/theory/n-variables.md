# Fluid states in n variables: a working note

Status: working note, October 2026. Parts 1 and 2 are proved below. Part 3 is a reduction whose
completeness still has to be shown, and Part 4 lists what is open. Nothing here is used by the
engine yet.

## Setting

Let P, Q be coprime real polynomials in x = (x1, ..., xn), n ≥ 2, with P(0) = Q(0) = 0, and
f = P/Q. Write S¹ = R ∪ {∞} for the circle (+∞ and −∞ are one point here; their signs are
recovered separately, as FMSL already does). The map

    Φ(x) = [P(x) : Q(x)] ∈ S¹

is defined and continuous off the indeterminacy set I = {P = Q = 0}; it equals f where Q ≠ 0 and
∞ where Q = 0 ≠ P. The fluid state of f at 0 is

    FS = ∩_{ε>0} closure Φ(B_ε \ ({0} ∪ I)).

(Points with Q = 0 ≠ P add only ∞, which is already a limit of f there, so this is the set of limit
values of f at 0 viewed on the circle.)

What prior work assumes, and this note does not: that Q vanishes only at 0 near 0. Without that
assumption, I can be a curve or a surface through 0 when n ≥ 3, and f is 0/0 at every point of it.
The 2016 Western Ontario work (Alvandi, Kazemi, Moreno Maza) calls this case work in progress;
the 2022 bivariate paper (Đinh, Guo, Nguyễn, Phạm) leaves non-isolated zeros to future work.

## Part 1. FS is a point, a closed arc, or the whole circle

P and Q coprime means I has codimension at least 2: over C their common zeros have dimension at
most n − 2, so the real set I has dimension at most n − 2. A closed semialgebraic set of
codimension at least 2 does not disconnect a connected open set, and B_ε \ {0} is connected for
n ≥ 2. So B_ε \ ({0} ∪ I) is connected, its image under Φ is connected, and so is the closure. FS
is a decreasing intersection of compact connected subsets of S¹, hence compact and connected.

## Part 2. A clean indeterminacy point forces the whole circle

**Theorem.** Suppose there are points p ∈ I \ {0} arbitrarily close to 0 at which dP(p) and dQ(p)
are linearly independent. Then FS = S¹ (every value, +∞ and −∞ included, is a limit value).

**Proof.** At such a point p, P and Q are part of a local coordinate system (implicit function
theorem): near p there are coordinates (u, v, w3, ..., wn) with P = u and Q = v. In every
neighbourhood of p, u/v takes every real value (take u = c v with v ≠ 0 small) and tends to ±∞
(take v → 0 with u fixed small, on either side). So every value of S¹, with either sign of ∞, is a
limit value of f at p. Choosing such p_k → 0 and, near each p_k, a point x_k within distance
|p_k|/k where f(x_k) is within 1/k of a given value c, gives x_k → 0 with f(x_k) → c. ∎

In two variables I \ {0} is empty near 0, so this case cannot arise; it is new to n ≥ 3, and it is
the typical case there: two random polynomials vanishing at 0 meet along a curve, usually with
independent differentials. When dP and dQ are dependent all along I (for example
f = x²/(x² + y²) in R³, where I is the z-axis and dP = 0 on it), the answer need not be the whole
circle: there FS = [0, 1].

## Part 3. Where the ends of a proper arc come from

Suppose FS is a proper arc. Pick a value c of S¹ outside FS and a Möbius map ψ sending c to ∞;
then g = ψ ∘ Φ is real-valued and bounded near 0 (off I), and FS corresponds to
[lim inf g, lim sup g]. On each small sphere S_r, the supremum of g over S_r \ I is either

  (a) attained at a point where the gradient of g is parallel to x, that is, on

        Γ = { x : rank [ x ; Q∇P − P∇Q ] ≤ 1 },

      which does not depend on c (ψ multiplies the gradient by a nonzero factor), and which also
      passes through the points where Q = 0 ≠ P, so infinity is treated like any other value; or

  (b) approached at points of I ∩ S_r, through values of f near I.

The same holds for the infimum. By the curve selection lemma, the contribution of (a) is the set of
limits of f along the half-branches of the curve part of Γ at 0 (plus finitely many critical
values on its higher-dimensional parts, where f is locally constant). The contribution of (b) is
the set J of limit values of f approached near I. So the ends of FS lie in

    E = { limits of f along half-branches of Γ } ∪ { critical values } ∪ J,

and FS is then decided by one test value in each gap between consecutive points of E: by Part 1 a
gap is either inside FS or outside it.

Γ itself is the discriminant variety of Alvandi, Kazemi and Moreno Maza, written for Φ instead of
f: using Q∇P − P∇Q rather than ∇(P/Q) is what lets it run through Q = 0.

## Part 4. What is open

1. **J, the limit values approached near I.** For a clean branch of I, Part 2 gives the whole
   circle. For a branch γ along which dP and dQ are dependent, f near a point p of γ is, in a
   slice transverse to γ, a 0/0 problem in n − 1 variables.

   *First guess, now refuted:* FS(p) is constant along γ near 0. It is not: in the exploration
   below (`benchmarks/part4_explore.py`, output in `docs/theory/part4_exploration.txt`), with γ the
   z-axis, FS at (0, 0, t) changes with t; for example [−∞, 0.00201] ∪ [0.24861, +∞] at t = 10⁻³
   and [−∞, −0.00199] ∪ [0.25136, +∞] at t = −10⁻³.

   *Refined conjecture:* the contribution of γ is the limit, as t → 0 from either side, of the
   fluid states FS(γ(t)) (a limit of closed sets on the circle), and the ends of FS(0) lie in that
   limit together with the limits along the branches of Γ. The slice fluid states FS(γ(t)) form a
   one-parameter family in n − 1 variables, which is what FMSL's moving-target engine computes in
   two variables (ends as exact functions of t), so the limit t → 0 is exact. This would give an
   induction on n: a non-clean branch in n variables becomes a moving target in n − 1.

   *Evidence (25 random functions, z-axis non-clean in each):* in all 7 cases whose answer is not
   simply the whole circle, FS(0) from the independent three-variable engine equals the limit of
   the slice fluid states, for example

   | f (abbreviated) | FS(0) | FS(0,0,t), t = 10⁻³, −10⁻³ |
   |---|---|---|
   | row 5 | [−∞, 0] ∪ [1/4, +∞] | [−∞, 0.00201] ∪ [0.24861, +∞]; [−∞, −0.00199] ∪ [0.25136, +∞] |
   | row 9 | {−∞} ∪ [1, +∞] | [0.999999, 999.626]; [−∞, −1000.38] ∪ [0.999999, +∞] |
   | row 16 | {−∞} ∪ [0, +∞] | [−∞, −1999] ∪ [0.000999937, +∞]; [−0.00100006, 2001] |

   (In row 9 the piece [−∞, −1000.38] for t < 0 shrinks to the single point −∞ as t → 0, which is
   why −∞ is an isolated value at the origin.) One case (row 10) has FS(0) = the whole circle while
   the slices along the z-axis give only {0}: there the whole circle comes from a different
   indeterminacy line, x = 0, y = z, which is clean, exactly as Part 2 predicts. These are 7
   cases, not a proof; they decide what to try to prove.

   *Stress test (`benchmarks/part4_test.py`, 60 functions, half built to break it):* this checks
   a stronger statement, that the tangential region alone (|(x, y)| = o(|z|)) adds nothing to
   the limit of the slice fluid states. Both sides are computed exactly. The tangential side is the
   fluid state at 0 of g(u, v, z) = f(zu, zv, z), and the slices are taken at t = ±10⁻⁶. Values
   above 10⁴ count as ∞, so slice ends that grow without bound as t → 0 are read as reaching ∞.

   | result | count |
   |---|---|
   | agree | 57 |
   | disagree | 1 |
   | engine unfinished | 2 |

   The disagreement is
   f = (−2xy(1+z) − 2zx³)/(−2x²z² − xy(2−z) − z³x). Its numerator and denominator share the
   factor x; the generator did not exclude that. After cancelling it,
   f = 2(y(1+z) + zx²)/(2xz² + y(2−z) + z³), and the z-axis is not in I at all. Each slice gives
   {0}, since f(0, 0, t) = 0/t³. Along x = 0, y = a z³, however, f → 2a/(2a+1), which takes every
   value as a varies (checked with exact arithmetic at z = 10⁻⁶). The denominator's zero surface,
   y ≈ −z³/2, hugs the axis at order z³, and the whole circle comes from there.

   So the stronger statement is false. Scales between the axis and the cone |(x, y)| ~ |z| can
   contribute, through branches of Γ (here through Q = 0 ≠ P) that are tangent to γ. The refined
   conjecture as stated above already includes the limits along the branches of Γ, and this example
   is consistent with it. What the example adds is that those branches must include the ones
   tangent to γ, which is the part to prove and to compute.
2. **Testing a gap value c.** c ∈ FS exactly when 0 is a limit value of
   (P − cQ)² / (P² + Q²), a function that is 0/0 only on I. An exact test is needed.
3. **Signs of ∞.** As in FMSL now: decided from the sign of Q on each side of its zero set near
   the branches where f is large.

## Inner bounds in any number of variables (`fmsl.witness`)

Curve selection says every limit value is reached along an analytic curve, and the leading
term of such a curve is c t^w. The module uses only curves of exactly that form,
x = (c1 t^w1, ..., cn t^wn). Each gives an exactly known limit whenever in_w Q(c) ≠ 0, and
segments in c give whole intervals. This is the one-term part of what `fmsl.space` does exactly
in three variables, without its zooming. So it proves values but can miss some.

In three variables it was checked against the exact engine on 50 functions: never outside, and
equal in 41. In four variables, 23,654 witnessed values were re-checked along their curves and
none was wrong (`benchmarks/witness_check.py`, `docs/theory/witness_check_output.txt`). The 9
three-variable misses all need curves with more than one term, such as x = y·sqrt(−z) on the zero
surface of x² + zy².

Version 0.32 adds one level of zooming, at rational zeros of the initial denominator. Of the 9,
6 are now found up to gaps narrower than 1e-11. The 3 left have a denominator vanishing on whole
coordinate planes (such as x²z). The repeated check: 47 of 50 equal to the exact answer, 8 of
those up to tiny gaps, and none outside. In four variables, 504,867 witnesses were re-checked and
none was wrong. What remains for an exact four-variable answer is an outer bound: a proof that no
other values occur. That is Part 3's reduction, made effective.

## Evidence so far

`benchmarks/part2_check.py`: 40 random pairs P = x a + y b, Q = x c + y d (a, b, c, d of degree at
most 1), which share the z-axis. Part 2's condition holds when det [[a, b], [c, d]] is not
identically 0 on the axis, which is decided exactly. FMSL's three-variable engine, which knows
nothing of this theorem, was run on each (60 s limit):

| | count |
|---|---|
| condition holds, engine proves the whole circle | 35 |
| condition holds, engine proves something else | 0 |
| condition holds, engine not proven within 60 s | 4 |
| condition fails (no prediction) | 1 |

No disagreement. The 4 cases the engine could not finish are answered by Part 2 at once, which is
the practical point: in three or more variables the theorem settles quickly what the zooming engine
finds hardest. This is a sanity check of a proved statement, not evidence for Parts 3 and 4; many
random three-variable functions have the whole circle as their answer anyway.
