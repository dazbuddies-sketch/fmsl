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
   circle. For a branch along which dP and dQ are dependent, f near a point p of the branch is,
   in a slice transverse to the branch, a 0/0 problem in n − 1 variables. The conjecture to check
   is that FS(p) is constant on each branch away from finitely many points, so that J is computed
   by a lower-dimensional fluid state, giving an induction on n.
2. **Testing a gap value c.** c ∈ FS exactly when 0 is a limit value of
   (P − cQ)² / (P² + Q²), a function that is 0/0 only on I. An exact test is needed.
3. **Signs of ∞.** As in FMSL now: decided from the sign of Q on each side of its zero set near
   the branches where f is large.

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
