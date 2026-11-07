# Celedon points hidden inside the zoom

*MJ Celedon. Module `fmsl.celedon`, used by `moving_profile(..., deep=True)`. Opt-in: the
default analysis is unchanged.*

## The problem

For a family f_t = N(x, y, t)/D(x, y, t), with coefficients polynomial in t, `fmsl.moving` finds
candidate Celedon points from the **top-level** Newton data: coefficients that vanish, and
discriminants and resultants of the initial forms. When the snapshot changes character only
after the engine zooms along a curve y = c₀|x|ᵃ + …, the change is invisible there. It was found
by sampling and reported as an approximate "change near …".

Example: f = x⁴ / ((y − x)² + (t+1)x⁴ − x·y³). Along y = x + c·x² the value tends to 1/(t + c²).
So the snapshot is [0, 1/t] for t > 0 and [−∞, 1/t] ∪ [0, +∞] for t < 0. The only Celedon point
is **t = 0**. The top-level data shows only t = −1, where the coefficient t + 1 of x⁴ vanishes but
nothing changes. The default analysis reports t = −1, a sampled "change near t ∈ (0, 10⁻⁶)", and,
because the interval (−1, ∞) is not split at 0, the end formula "v = 1" (true only at t = 1).
Deep mode reports exactly t = 0, no unexplained change, and the end formula t·v − 1 = 0
(v = 1/t) on both sides.

## The method: the zoom with t kept as a symbol

At a tie weight a = p/q and a side σ = ±1 of x, let f_N(c, t), f_D(c, t) be the initial forms
(polynomials in c over Q[t]). The engine zooms at the real roots c₀ that are repeated in f_N or
f_D, or shared by both. These are the roots of

  G = gcd_c(f_N, ∂f_N),  gcd_c(f_D, ∂f_D),  gcd_c(f_N, f_D)   (gcds over Q(t)).

1. **Generic validity.** Let h be the squarefree part (over Q(t)) of f_N·f_D with the factor c
   removed. Add to the candidates the real roots of Res_c(h, ∂h), of h's leading coefficient in c,
   and of its constant term.
2. **Zoom roots in Q(t).** For each G: pick a rational t₀ where G(·, t₀) keeps its degree and is
   squarefree; at each rational root r₀ of G(·, t₀), expand the root as a power series in t − t₀
   (Newton lifting), take Padé approximants of degree ≤ deg_t G, and accept c₀ = P/Q only if
   G(P/Q, t) = 0 **identically** (exact polynomial check). If fewer than deg_c G roots are found
   this way, the branch is reported as **unsupported**.
3. **Substitution.** x = σ sᵠ, y = sᵖ (P(t) + w)/Q(t), and both N and D multiplied by Q^K (same
   K, so the ratio is unchanged). This gives a polynomial family N₁, D₁ in (s, w, t). Add the
   real roots of Q, and the top-level candidates of (N₁, D₁) (`candidate_times`).
4. **Recurse** on (N₁, D₁) to a fixed depth (default 3). Branches beyond it are reported as
   unsupported.

5. **Arcs of different branches.** Collect every end relation E(t, v) = 0, meaning the boundary
   and critical values at every node of the top level and of every zoomed family. Add the real
   roots of three kinds of polynomial in t:
   * Res_v(E_i, E_j) for every pair, after dividing out a shared factor;
   * the discriminant in v of each E_i;
   * the leading coefficient in v of each E_i.

   These are the times where two ends meet, so arcs can merge or separate, or where an end escapes
   to infinity. If some end of a reported interval has no exact relation, deep mode says so in
   `deep_unsupported`.

In deep mode, end formulas are also taken from the zoomed families. In **both** modes, an end
formula is now accepted only if the matching end also satisfies it at the other sample times of
the interval. That removes formulas that fit at a single time by coincidence.

In the default mode, a family that zooms along a curve for generic t now carries a **note** saying
that changes hidden in the zoom are not searched there, and can show up displaced to the nearest
visible candidate. It changes no result. In the stress test, the note appeared on all 39 families
where the default mode made an error, and on 15 of the 98 families where it made none.

## Statement

**Proposition.** Let 𝒞 be the set of real times returned by `all_candidate_times(N, D, depth)`, and
suppose its list of unsupported branches is empty. Let I be an open interval that contains no
point of 𝒞. Then for t ∈ I, along every zoom branch of depth ≤ `depth`:

(a) the set of tie weights, and the degrees of the initial forms in c, are the same for all t ∈ I;

(b) the real roots of f_N·f_D (other than c = 0) form a fixed number of continuous functions of t,
which never meet, never leave the real line, and keep their multiplicities in f_N and in f_D; in
particular the zoom points are the same branches c₀(t) = P(t)/Q(t) for all t ∈ I, with Q ≠ 0;

(c) the critical points of the ratio f_N/f_D never collide with each other or with a pole.

So the engine's zoom tree has the same shape for every t ∈ I, and at every node the local data it
reads (which ends are infinite, the order of critical values and poles, where the numerator
changes sign) has the same combinatorial type.

**Proof.** (a) A tie weight or a degree in c can change only where some coefficient, a polynomial
in t, vanishes. At the top level these are the coefficient roots in `candidate_times`. Inside a
zoom they are the coefficient roots of (N₁, D₁), which step 3 adds. (b) Let h be the squarefree
part of f_N·f_D over Q(t). On an interval where h's leading coefficient and constant term in c
do not vanish, its degree is fixed and no root goes to 0 or ∞. Where Res_c(h, ∂h) does not vanish,
h(·, t) is squarefree, so its roots are distinct and depend continuously on t (by the implicit
function theorem). A real root can become complex only by meeting its conjugate, which is a
double root, excluded. The multiplicity of a root in f_N or f_D can jump only if two distinct
roots of h merge, also excluded. The gcd G specializes correctly at every t where no two roots
of h merge. So the zoom roots are exactly the branches found in step 2. The unsupported list is
empty, so every one of them is a rational function P/Q, and Q has no root in I by step 3. (c) These
are the "critical points collide" and "critical point meets a pole" resultants of
`candidate_times`, applied at the top level and to every (N₁, D₁). ∎

**Theorem (whole snapshot).** Let 𝒞 now include the arc candidates of step 5. Suppose
`deep_unsupported` is empty, so all zoom roots are rational in t, the depth suffices, and every
end has an exact relation. Then on every open interval I that contains no point of 𝒞, the snapshot
has the same type for all t ∈ I. Hence every Celedon point lies in 𝒞, and the Celedon points
`moving_profile(deep=True)` reports are exactly the points of 𝒞 where the type actually changes.

**Proof.** The snapshot is the union of the arcs produced at the nodes of the zoom tree. By the
Proposition, the tree, and the set of arcs each node produces, keep the same combinatorial
description on I. Each finite end of an arc is a real root, in v, of one of the relations E_i(t, ·),
and moves continuously with t. By step 5, on I no end becomes infinite (leading coefficients),
no two ends coming from the same relation meet (discriminants), and no two ends coming from
different relations meet (resultants). So the order of all ends on the extended line, with ±∞
included, is the same for every t ∈ I. The union of a fixed family of arcs depends only on the
order of their ends. So the type of the union, meaning which arcs overlap, merge or stay apart,
is constant on I. ∎

Example: f = x⁴/((y − x)² − x⁴) + t·x²/(x² + y²). The ordinary directions give [0, t] for t > 0,
and [t, 0] for t < 0. The curve near y = x gives [−∞, t/2 − 1] ∪ (t/2, +∞]. The pieces merge into
the whole line at **t = 2**, where t/2 − 1 meets 0, and at **t = −2**, where t/2 − 1 meets t. At
t = −2 nothing changes inside any single branch, so before step 5 this point was only found
approximately. Now both are exact. The default mode still gives t = −2 only approximately, and
reports a displaced point at t = 1 with the note.

## Zoom curves that are irrational in t (`fmsl.celedon_alg`)

When a zoom root is not a rational function of t, let h(u, t) be the factor of G that is left after
the roots in Q(t) are removed. The substitution x = σ·s^q, y = s^p·(u + w) gives a family whose
coefficients are polynomials in (t, u), to be read on the curve h(u, t) = 0. The method follows
the same steps as `candidate_times`, with three changes:

1. **Support on the branch.** A coefficient is dropped when h divides it (gcd over Q(t)). If it
   vanishes on only one factor of h, h is split and each part is searched separately.
2. **Resultants in c with coefficients in Q[t, u].** These are computed exactly: at enough values
   of u where the leading coefficients survive, then interpolated in u and checked at one more value.
3. **Elimination of u.** Each condition P(t, u) becomes R(t) = Res_u(P, h). Every t where P vanishes
   at some root u of h(·, t) is a root of R, so nothing on the branch is missed. Times from conjugate
   branches may be added, which is harmless.

The branch's end relations, meaning boundary values and critical values with u eliminated the
same way, are added to the arc-merging step and to the end formulas.

Example: x⁸/((y² − t·x²)² + (t − 1 + t²)x⁸ − x⁴y⁴). Along y = √t·x the term x⁴y⁴ equals t²x⁸ and
cancels the hidden part of the coefficient, so the value tends to 1/(4t·c² + t − 1). The snapshot is
bounded for t > 1 and the whole line for 0 < t < 1. The Celedon point **t = 1** is invisible at the
top level, where the candidates are 0 and the roots of t² + t − 1. Deep mode finds it exactly. The
default mode reports a displaced point at 0.618, the echo, and leaves t = 1 approximate. The same
works for a constant irrational curve (y = √2·x) and for y = √(t+1)·x.

## What is not covered (stated plainly)

* **A further zoom inside an irrational branch.** Only one irrational zoom level is searched.
  Deeper zooms inside it are listed in `deep_unsupported` and fall back to sampling.
* **Zoom trees deeper than the depth limit** are also listed as unsupported.
* **Ends without an exact relation** are listed as unsupported; merges involving them are checked
  by sampling only.

## Relation to known results

That a polynomial family has only finitely many such change points follows from the
Tarski–Seidenberg theorem (the snapshot is first-order definable over the reals), and general
quantifier elimination or cylindrical algebraic decomposition can in principle compute them, at
doubly exponential cost. What is offered here is an exact, local method that follows the engine's
own zoom, with a proof of what it covers and an explicit list of what it does not. I did not find
prior work computing exact limit sets at a 0/0 point for families with a parameter. The search
was brief, so this is provisional.

## End relations at a level with a shared root (both modes)

When f_N and f_D share a factor, which happens along every curve the engine zooms on, the
critical-value resultant Res_c(g, v·f_D − f_N) vanishes identically. The end relations of that
level were then silently lost. Without them, step 5 could not see arcs merge, and deep mode
reported 5/16 where the true merge is at 5/18. The common factor is now cancelled first. The ratio
is unchanged, so this only adds correct relations.

## Ends that switch branch (both modes)

The shape of the snapshot can stay the same while one of its ends passes from one branch to
another. For example, x⁴/((y + x)² + x⁴) + t·y²/(x² + y²) is a single arc for every t, but its right
end is t/2 + 1 for 0 < t < 2 and t for t > 2. Such a time is not a Celedon point. Stretches on
either side used to be joined and given a single end formula, which was wrong on part of it. Now
they stay separate, and the boundary is marked "end switches branch" (`end_switch_at`). Deep mode
finds every such switch exactly, since step 5 includes it. Default mode finds only those among its
own candidates.

## Unproven snapshots (both modes)

At an irrational time where the zoom curve is irrational too (y = √2·x), the two-variable engine
returns a partial snapshot marked `exact = False`. The moving layer used to treat it as a fact. That
produced false Celedon points (for example at t = ∛(1/2) for
(x⁴ − x·y³)/((y − t·x)² + (t + 1 + 2t³)x⁴ − 2x·y³), where the true snapshot is [0, 0.2788], not {0}).
Now an unproven snapshot is never on its own evidence of a change, and its text is labeled
"(partial)". See `benchmarks/results/celedon_stress/ERROR_LOG.md`.

## Checks

* `tests/test_celedon_deep.py`: the hidden example (exact t = 0, end formula t·v = 1), the default
  mode unchanged apart from the rejected coincidental end formula, roots in Q(t) (t, 2, 1/t),
  irrational branches reported (c² + t, c² − t), the substitution preserving the ratio, and a zoom
  root that moves with t (y = t·x, hidden point t = 2).
* `benchmarks/celedon_deep_check.py`: classic, random, and hidden-zoom families, in both modes,
  with dense sampling (12 times per interval) of the deep result.
* The full existing test suite passes unchanged.
