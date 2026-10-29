"""Euclidean distances in two and three variables (fmsl.distance, through fmsl.verify)."""
import unittest

from fmsl.distance import reduce_distance, DistanceError
from fmsl.verify import limit_at_origin, limit_at, check_claim

R3 = "sqrt(x^2+y^2+z^2)"


class TestDistance(unittest.TestCase):
    def test_two_variables(self):
        for e, s in [("a*b/sqrt(a^2+b^2)", "{0}"), ("(sqrt(x^2+y^2)-x)/sqrt(x^2+y^2)", "[0, 2]"),
                     ("(sqrt(x^2+y^2)-x)/y", "[-oo, +oo]")]:
            with self.subTest(e=e):
                r = limit_at_origin(e); self.assertTrue(r["proven"]); self.assertEqual(r["set"], s)

    def test_three_variables(self):
        for e, s in [(f"x/{R3}", "[-1, 1]"), (f"(x^2+y^2)/{R3}", "{0}"), (f"1/{R3}", "{+oo}"),
                     (f"({R3}-x)/{R3}", "[0, 2]"), ("z/sqrt(x^2+y^2)", "[-oo, +oo]")]:
            with self.subTest(e=e):
                r = limit_at_origin(e); self.assertTrue(r["proven"]); self.assertEqual(r["set"], s)

    def test_algebraic_ends(self):
        (lo, hi), = limit_at_origin(f"x*y*z/{R3}^3")["components"]   # +-1/(3 sqrt 3)
        self.assertAlmostEqual(hi, 1 / 27 ** 0.5, places=12); self.assertAlmostEqual(lo, -1 / 27 ** 0.5, places=12)

    def test_other_points(self):
        self.assertEqual(limit_at("(x-1)*y/sqrt((x-1)^2+y^2)", [1, 0])["set"], "{0}")
        self.assertEqual(limit_at("(x-1)/sqrt(x^2-2*x+1+(y-2)^2+z^2)", [1, 2, 0])["set"], "[-1, 1]")

    def test_distance_becomes_a_polynomial(self):
        # under the map the distance itself is the polynomial u^2 + v^2 + w^2
        g = reduce_distance(f"{R3}", ["x", "y", "z"])
        self.assertEqual(g, "(x^2+y^2+z^2)")

    def test_refusals(self):
        for e in ["x/sqrt(x^2+2*y^2)", "sqrt(x^2+y^2)/sqrt(y^2+z^2)"]:
            with self.subTest(e=e), self.assertRaises(DistanceError):
                limit_at_origin(e)

    def test_check_claim(self):
        self.assertTrue(check_claim(f"x*y*z/{R3}^3", "[-0.19245, 0.19245]")["correct"])
        self.assertFalse(check_claim(f"x/{R3}", "0")["correct"])


if __name__ == "__main__":
    unittest.main()
