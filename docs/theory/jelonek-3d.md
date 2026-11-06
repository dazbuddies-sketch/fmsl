# Jelonek sets in three variables: a working note (November 2026)

Question (Z. Jelonek, November 2026): find polynomial maps F: R^3 -> R^3 whose Jelonek set
(non-properness set) S_F has codimension 2, preferably with a nowhere-vanishing Jacobian.
For a local diffeomorphism, a nonempty S_F of codimension 2 would contradict Jelonek's real
Jacobian conjecture, which is open exactly in codimension 2 for n >= 3 (it is known for
codimension >= 3 and for n = 2; Pinchuk's map shows codimension 1 is not enough).

## Tool: `fmsl.jelonek_n`

- Inner bound: points of S_F, each proven by an explicit arc to infinity (monomial arcs plus one
  level of zooming at rational points).
- Outer bound: for each coordinate, an exact relation R(x_i, F(x)) = 0 found by linear algebra
  modulo 2^61 - 1, reconstructed over Q and checked exactly. Its leading coefficient a(u)
  vanishes at every point of S_F reached with x_i unbounded, so S_F is contained in the union of
  the V(a_i), in each of two coordinate systems.
- Dimension: random lines are checked exactly (Sturm) against the outer bound. A surface is met
  by a positive share of lines, a curve almost never. This is evidence, not proof.

## A simple observation: triangular maps cannot answer the question

For F = (x, y, h(x, y, z)), det J_F = h_z. If h_z >= 0 everywhere and vanishes on no segment of
a vertical line, then h is strictly increasing on every vertical line, so it is a polynomial of
odd degree in z and maps that line onto R. Then F is a continuous bijection of R^3, hence a
homeomorphism by invariance of domain, hence proper: S_F is empty. So a triangular map with nonempty S_F has h_z
changing sign, and det J_F vanishes on a surface. The same holds after composing with polynomial
automorphisms.

## Search (`benchmarks/jelonek3_search.py`, output `jelonek3_search_output.txt`)

300 maps of five kinds:
- triangular (x, y, A(x, y) z^d + ...) with A vanishing only at 0;
- (x, g, h) with g and h random;
- quadratic perturbations of the identity;
- random cubic maps;
- triangular codimension-2 maps composed with polynomial automorphisms (a test of the tool on
  maps that are not visibly triangular).

| result | count |
|---|---|
| codimension-2 candidates (proven points; no line of 200 meets the outer bound) | 80 |
| ...of which det J vanishes on a surface | 80 |
| S_F has proven points; outer bound contains a surface | 129 |
| no proven point, outer bound has dimension <= 1 by the line test | 59 |
| no proven point, outer bound contains a surface | 3 |
| not finished within 30 s | 29 |

All 80 candidates are triangular (45) or triangular in disguise (35). Every one has
S_F = the z-axis line (or the image of it). No random map of the other three kinds gave a
candidate. So the search found many codimension-2 Jelonek sets, as expected, and nothing with
the preferred Jacobian condition, which is consistent with the conjecture. A failed search
proves nothing.

## The 2026 counterexample to the Jacobian conjecture

F = ((1+xy)^3 z + y^2(1+xy)(4+3xy), y + 3x(1+xy)^2 z + 3xy^2(4+3xy), 2x - 3x^2 y - x^3 z),
with det J_F = -2 (the counterexample of July 2026, "A counterexample to the Jacobian
Conjecture", ulam.ai).

- Outer bound, exact. The relations for y and z have constant leading coefficients, so along any
  sequence to infinity with F bounded, y and z stay bounded and only x grows. The relation for x
  has leading coefficient (times 27)
        A(u) = 27 u1^2 u3^2 - 18 u1 u2 u3 + u2^3 u3 + 16 u1 - u2^2,
  so the real Jelonek set lies in the surface A = 0.
- Inner bound, exact, two-dimensional. For any real a, b, g, h, along the arc
        x = 1/s,  y = a + b s + g s^2,  z = s(-3a + (2 - 3b) s + h s^2),   s -> 0+,
  F is a polynomial in s plus the terms below, so F tends to
        L(a, lambda) = (a^2 - a^3 lambda,  4a - 3a^2 lambda,  lambda),   lambda = -(h + 3g).
  The derivative of L at (a, lambda) = (1, 0) has rank 2, so these limits fill a two-dimensional
  piece of the Jelonek set, and A(L(a, lambda)) = 0 identically. Checked three ways: the Laurent
  expansion (the s^-1 coefficients of all three components are multiples of e + 3b - 2), exact
  rational evaluation along 300 random arcs at s = 10^-30 (error below 10^-24, scaling like s),
  and symbolic expansion of A(L) = 0.
- The automatic search (three levels of zooming) separately proves 901 points with arcs, 542 off
  the u3-axis, all on A = 0; each of its 14,184 arcs was re-checked exactly as a Laurent
  polynomial by separate code. Each individual arc family it uses is only one-dimensional, so the
  two-dimensional statement rests on the explicit family above.
- The u3-axis lies both in the Jelonek set and in the image: F(c/2, 0, 0) = (0, 0, c). (An
  earlier version of this note, and an email of 6 November 2026, wrongly said that (0, 0, c),
  c != 0, is not in the image; that came from an algebra slip and has been corrected.)

## Second search: critical values (November 6)

Z. Jelonek (6 November 2026): a map with S_F of codimension 2 whose critical values miss S_F
would also be interesting. 400 more maps (`jelonek3_search2_output.txt`), adding three kinds that
are not visibly triangular. Results: 72 codimension-2 candidates, every one with S_F a line and
det J vanishing on a surface, and every one triangular up to polynomial automorphisms of the
source or the target. The two genuinely non-triangular kinds gave none.

Critical values. The first test in the search script (Newton's method on F1 = F2 = det J = 0)
was unreliable: it misses solutions that form a curve, and reported "no solution" for all 72.
It was replaced by `benchmarks/jelonek3_critical_check.py`, checked two ways:
- Levenberg-Marquardt minimisation of F1^2 + F2^2 + det^2 (handles curves of solutions);
- for plainly triangular maps (x, y, h), an exact criterion: the critical values meet the line
  u1 = u2 = 0 exactly when h_z(0, 0, z) has a real root or vanishes identically (Sturm's
  theorem).
The two agree on all 42 plainly triangular maps. Of the 72: 31 have critical values meeting
S_F and 41 do not (like Example 3.4 of Braun-Jelonek-Michalska). Nothing found goes beyond
triangular maps up to automorphisms.
