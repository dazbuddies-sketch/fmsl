import unittest
from fractions import Fraction as Fr

from fmsl.moving import moving_profile, parse_moving, snapshot, endpoint_formulas, _solve_text


class MovingTargets(unittest.TestCase):
    def test_special_time_and_intervals(self):
        r = moving_profile("(x^2-t*y^2)/(x^2+y^2)")
        self.assertEqual([s["t"] for s in r["special"]], [-1])
        self.assertEqual(r["special"][0]["snapshot"], "{1}")
        self.assertEqual(r["special"][0]["moving"], "{1}")
        self.assertEqual(r["unexplained"], [])
        texts = [iv["snapshot"] for iv in r["intervals"]]
        self.assertIn("[1, 3]", texts)                       # t = -3

    def test_endpoint_formulas(self):
        N, D = parse_moving("x^2*y/(x^4+t*y^2)")
        s = snapshot(N, D, Fr(2))
        forms = {_solve_text(E) for _, E in endpoint_formulas(N, D, Fr(2), s["comps"]) if E}
        self.assertEqual(forms, {"4*t*v^2 - 1 = 0  (v = the end)"})   # ends +-1/(2 sqrt t)
        N, D = parse_moving("(x^2-t*y^2)/(x^2+y^2)")
        s = snapshot(N, D, Fr(2))
        forms = {_solve_text(E) for _, E in endpoint_formulas(N, D, Fr(2), s["comps"]) if E}
        self.assertEqual(forms, {"v = -t", "v = 1"})

    def test_irrational_special_time(self):
        r = moving_profile("(x*y+t*x^2)/(x^2+(t-1)*y^2)", joint=False)
        ts = [float(s["t"]) for s in r["special"]]
        self.assertAlmostEqual(ts[0], -0.7548776662466927, places=12)   # root of t^3 - t^2 + 1
        self.assertEqual(r["special"][0]["snapshot"], "[-oo, +oo]")

    def test_moving_path(self):
        # the target slides along (t, t^2): a 0/0 point only at t = 0
        r = moving_profile("x^2*y/(x^4+y^2)", path=("t", "t^2"))
        self.assertEqual([s["t"] for s in r["special"]], [0])
        self.assertEqual(r["special"][0]["moving"], "[-1/2, 1/2]")


if __name__ == "__main__":
    unittest.main()
