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

In deep mode, end formulas are also taken from the zoomed families. In **both** modes, an end
formula is now accepted only if the matching end also satisfies it at the other sample times of
the interval. That removes formulas that fit at a single time by coincidence.

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

## What is not covered (stated plainly)

* **Irrational zoom roots** (for example c₀ = √t). These branches are listed in
  `deep_unsupported` and fall back to sampling, exactly as before.
* **Zoom trees deeper than the depth limit** are also listed as unsupported.
* **Ends from different branches meeting.** Two arcs coming from different zoom branches can
  merge when their ends cross (for example 1/t on one branch meets 2 on another, at t = 1/2).
  The Proposition fixes each branch's local type, not the way the arcs of different branches
  overlap. The same gap exists at the top level. A change of this kind is still found by sampling
  and reported as "change near …". Closing it needs resultants between the end relations
  E(t, v) of different branches. That is the next step.

## Relation to known results

That a polynomial family has only finitely many such change points follows from the
Tarski–Seidenberg theorem (the snapshot is first-order definable over the reals), and general
quantifier elimination or cylindrical algebraic decomposition can in principle compute them, at
doubly exponential cost. What is offered here is an exact, local method that follows the engine's
own zoom, with a proof of what it covers and an explicit list of what it does not. I did not find
prior work computing exact limit sets at a 0/0 point for families with a parameter. The search
was brief, so this is provisional.

## Checks

* `tests/test_celedon_deep.py`: the hidden example (exact t = 0, end formula t·v = 1), the default
  mode unchanged apart from the rejected coincidental end formula, roots in Q(t) (t, 2, 1/t),
  irrational branches reported (c² + t, c² − t), the substitution preserving the ratio, and a zoom
  root that moves with t (y = t·x, hidden point t = 2).
* `benchmarks/celedon_deep_check.py`: classic, random, and hidden-zoom families, in both modes,
  with dense sampling (12 times per interval) of the deep result.
* The full existing test suite passes unchanged.
