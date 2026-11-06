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
