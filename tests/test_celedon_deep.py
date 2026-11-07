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
        t1 = r["intervals"][1]["sample"]
        self.assertEqual([iv["snapshot"] for iv in r["intervals"]],
                         ["[-oo, -1/3] u [0, +oo]", f"[0, {1 / t1}]"])     # [0, 1/t] for t > 0
        # the moving end is v = 1/t, i.e. t*v - 1 = 0, on both sides
        tv = {(1, 1): Fr(1), (0, 0): Fr(-1)}
        self.assertEqual([E for _, E in r["intervals"][0]["ends"]][0], tv)
        self.assertEqual([E for _, E in r["intervals"][1]["ends"]][1], tv)

    def test_default_mode_unchanged_points(self):
        # default mode still reports only what the top-level data shows (the hidden point is
        # left as an approximate change), but no longer prints a coincidental end formula
        r = moving_profile(HIDDEN, deep=False)
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
        r = moving_profile(f, joint=False, deep=False)
        self.assertNotIn(0.7937, [round(float(s["t"]), 4) for s in r["celedon_points"]])

    def test_arcs_from_different_branches_merge(self):
        # ordinary directions give [0, t] (t > 0) or [t, 0] (t < 0); the hidden curve near y = x
        # gives t/2 + 1/(c^2 - 1), i.e. [-oo, t/2 - 1] u (t/2, +oo].  The two pieces merge into the
        # whole line when t/2 - 1 reaches 0 (t = 2) or reaches t (t = -2).  Nothing changes inside
        # a single branch at t = -2: only the way the arcs combine.
        f = "x^4/((y-x)^2-x^4) + t*x^2/(x^2+y^2)"
        r = moving_profile(f, deep=True, joint=False)
        self.assertEqual([s["t"] for s in r["celedon_points"]], [-2, 2])
        self.assertEqual(r["unexplained"], [])
        self.assertTrue(any("ends of two branches meet" in w
                            for s in r["celedon_points"] if s["t"] == -2 for w in s["why"]))
        # the default mode still sees t = -2 only approximately (unchanged behaviour)
        r0 = moving_profile(f, joint=False, deep=False)
        self.assertEqual(len(r0["unexplained"]), 1)

    def test_merge_needs_ends_of_a_level_with_a_shared_root(self):
        # ordinary arc [-|t|, |t|] (2t*c/(1+c^2)); hidden curve y ~ 2x gives 4t/5 + 1/(c^2-2).
        # The arcs merge at t = 5/18 and t = -5/2.  The ends +-t used to be lost because the
        # level's initial forms share the factor (c-2)^2, which killed the end-formula resultant.
        r = moving_profile("x^4/((y-(2)*x)^2-2*x^4) + 2*t*x*y/(x^2+y^2)", deep=True, joint=False)
        self.assertEqual([s["t"] for s in r["celedon_points"]], [Fr(-5, 2), Fr(5, 18)])
        self.assertEqual(r["unexplained"], [])
        self.assertEqual(r["deep_unsupported"], [])

    def test_end_switch_is_not_a_celedon_point(self):
        # one arc for every t: [0,t] u [t/2, t/2+1] (t>0), [t,0] u [t/2, t/2+1] (t<0); the ends
        # switch branch at t = -2, 0, 2 without a change of shape
        r = moving_profile("x^4/((y-(-1)*x)^2+1*x^4) + t*y^2/(x^2+y^2)", deep=True, joint=False)
        self.assertEqual(r["celedon_points"], [])
        self.assertEqual([iv.get("end_switch_at") for iv in r["intervals"]], [None, -2, 0, 2])
        right = {(0, 1): Fr(2), (0, 0): Fr(-2), (1, 0): Fr(-1)}           # 2v - 2 - t = 0
        self.assertEqual(r["intervals"][1]["ends"][1][1], right)
        self.assertEqual(r["intervals"][3]["ends"][1][1], {(0, 1): Fr(-1), (1, 0): Fr(1)})   # v = t

    def test_hidden_on_irrational_branches(self):
        # along y = u x (u^2 = p), x^4 y^4 = p^2 x^8 cancels the hidden part of the x^8
        # coefficient, leaving (t - r): the value is 1/(4 p c^2 + t - r) on the branch
        cases = [("x^8/((y^2-t*x^2)^2 + (t-1+t^2)*x^8 - x^4*y^4)", [0, 1]),            # u = sqrt(t)
                 ("x^8/((y^2-(t+1)*x^2)^2 + (t-2+(t+1)^2)*x^8 - x^4*y^4)", [-1, 2]),    # u = sqrt(t+1)
                 ("(x^8+t*x^4*y^4)/((y^2-2*x^2)^2 + (t-1+4)*x^8 - x^4*y^4)", [Fr(-1, 4), 1])]  # u = sqrt 2
        for f, expected in cases:
            r = moving_profile(f, deep=True, joint=False)
            self.assertEqual([s["t"] for s in r["celedon_points"]], expected, f)
            self.assertEqual(r["unexplained"], [], f)
            self.assertEqual(r["deep_unsupported"], [], f)

    def test_end_zero_has_a_formula(self):
        # after the common factor x cancels, f tends to 0 along every path; the end v = 0 is never
        # a boundary or critical value, and used to be reported as an end without a formula
        f = "((-2)*(t-1)*x^3*y^2)/((-1)*(t-2)*x^3*y^0 + (-2)*(t-1)*x^4*y^4)"
        r = moving_profile(f, joint=False)
        self.assertEqual([s["t"] for s in r["celedon_points"]], [2])
        self.assertEqual(r["deep_unsupported"], [])
        self.assertTrue(any(E == {(0, 1): Fr(1)} for iv in r["intervals"] for _, E in iv["ends"]))

    def test_deep_is_the_default(self):
        r = moving_profile(HIDDEN)
        self.assertIn("deep_unsupported", r)
        self.assertEqual([s["t"] for s in r["celedon_points"]], [0])

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
