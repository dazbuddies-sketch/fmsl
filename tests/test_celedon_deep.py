import unittest
from fractions import Fraction as Fr

from fmsl.moving import moving_profile, parse_moving
from fmsl.celedon import rational_roots_in_t, zoom_substitute, deep_candidate_times

HIDDEN = "x^4/((y-x)^2+(t+1)*x^4-x*y^3)"      # along y = x + c x^2 the value is 1/(t + c^2)


def _ev(F, a, b, t):
    return sum(v * a ** i * b ** j * t ** k for (i, j, k), v in F.items())


class DeepCeledon(unittest.TestCase):
    def test_hidden_point_found_exactly(self):
        r = moving_profile(HIDDEN, deep=True)
        self.assertEqual([s["t"] for s in r["celedon_points"]], [0])
        self.assertEqual(r["unexplained"], [])
        self.assertEqual(r["deep_unsupported"], [])
        self.assertEqual([iv["snapshot"] for iv in r["intervals"]],
                         ["[-oo, -1/3] u [0, +oo]", "[0, 1/2]"])
        # the moving end is v = 1/t, i.e. t*v - 1 = 0, on both sides
        tv = {(1, 1): Fr(1), (0, 0): Fr(-1)}
        self.assertEqual([E for _, E in r["intervals"][0]["ends"]][0], tv)
        self.assertEqual([E for _, E in r["intervals"][1]["ends"]][1], tv)

    def test_default_mode_unchanged_points(self):
        # default mode still reports only what the top-level data shows (the hidden point is
        # left as an approximate change), but no longer prints a coincidental end formula
        r = moving_profile(HIDDEN)
        self.assertEqual([s["t"] for s in r["celedon_points"]], [-1])
        self.assertEqual(len(r["unexplained"]), 1)
        ends = r["intervals"][1]["ends"]
        self.assertEqual(ends[1][1], None)          # was "v = 1", which is wrong for t != 1

    def test_rational_roots_in_t(self):
        # (c - t)^2 (c - 2): roots t and 2, both rational in t
        from fmsl.moving import _bmul
        from fmsl.algebraic import pmul, ptrim
        g = _bmul(_bmul({(0, 1): Fr(1), (1, 0): Fr(-1)}, {(0, 1): Fr(1), (1, 0): Fr(-1)}),
                  {(0, 1): Fr(1), (0, 0): Fr(-2)})
        roots, complete = rational_roots_in_t(g)
        self.assertTrue(complete)
        self.assertEqual(sorted((tuple(P), tuple(Q)) for P, Q in roots),
                         sorted([((Fr(0), Fr(1)), (Fr(1),)), ((Fr(2),), (Fr(1),))]))
        # (t c - 1): root 1/t
        roots, complete = rational_roots_in_t({(1, 1): Fr(1), (0, 0): Fr(-1)})
        self.assertTrue(complete)
        P, Q = roots[0]
        self.assertEqual(ptrim(pmul(P, [Fr(0), Fr(1)])), ptrim(Q))        # P/Q = 1/t

    def test_irrational_branch_is_reported(self):
        # c^2 + t: real roots +-sqrt(-t) for t < 0, complex at t = 1 -- must NOT count as covered
        _, complete = rational_roots_in_t({(0, 2): Fr(1), (1, 0): Fr(1)})
        self.assertFalse(complete)
        _, complete = rational_roots_in_t({(0, 2): Fr(1), (1, 0): Fr(-1)})   # c^2 - t
        self.assertFalse(complete)

    def test_substitution_keeps_the_ratio(self):
        N, D = parse_moving(HIDDEN)
        P, Q = [Fr(1)], [Fr(1)]
        N1, D1 = zoom_substitute(N, D, Fr(1), 1, P, Q)
        for s, w, t in [(Fr(1, 3), Fr(1, 5), Fr(2)), (Fr(1, 7), Fr(-1, 2), Fr(-3))]:
            x, y = s, s * (1 + w)
            self.assertEqual(_ev(N1, s, w, t) / _ev(D1, s, w, t), _ev(N, x, y, t) / _ev(D, x, y, t))

    def test_unproven_snapshot_is_not_a_celedon_point(self):
        # at the algebraic time t0 = cbrt(1/2) ~ 0.7937 the zoom curve y = t0*x is irrational, the
        # two-variable engine returns a PARTIAL snapshot ({0}, exact=False), while the true snapshot
        # is [0, (1-t0^3)/(t0+1)] ~ [0, 0.2788], the same type as on both sides.  It used to be
        # reported as a Celedon point.
        f = "(x^4-x*y^3)/((y-(t)*x)^2 + (t-(-1)+2*(t)^3)*x^4 - 2*x*y^3)"
        r = moving_profile(f, deep=True, joint=False)
        self.assertEqual([s["t"] for s in r["celedon_points"]], [-1, 1])
        self.assertEqual(r["unexplained"], [])
        r = moving_profile(f, joint=False)
        self.assertNotIn(0.7937, [round(float(s["t"]), 4) for s in r["celedon_points"]])

    def test_t_dependent_zoom_root(self):
        # zoom root y = t*x; hidden coefficient (t - 2) at the second level
        f = "x^4/((y-t*x)^2+(t-2+t^3)*x^4-x*y^3)"
        N, D = parse_moving(f)
        cands, unsupported = deep_candidate_times(N, D)
        self.assertIn(Fr(2), [t for t, _ in cands])
        r = moving_profile(f, deep=True, joint=False)
        self.assertIn(Fr(2), [s["t"] for s in r["celedon_points"]])
        self.assertEqual(r["unexplained"], [])


if __name__ == "__main__":
    unittest.main()
