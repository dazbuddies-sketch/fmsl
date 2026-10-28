"""Limits with the distance r = sqrt(x^2 + y^2) (fmsl.radial)."""
import unittest

from fmsl.radial import radial_limit
from fmsl.verify import limit_at_origin


class TestRadial(unittest.TestCase):
    def test_sets(self):
        for text, expected in [("x*y/r", "{0}"), ("x*y/sqrt(x^2+y^2)", "{0}"), ("(r-x)/y", "[-oo, +oo]"),
                               ("x^2/(r^2+r^3)", "[0, 1]"), ("(r-x)/r", "[0, 2]"), ("1/r", "{+oo}")]:
            with self.subTest(text=text):
                r = radial_limit(text)
                self.assertTrue(r["proven"])
                self.assertEqual(r["set"], expected)

    def test_limit_value(self):
        r = radial_limit("(r^2-x^2)/y^2")
        self.assertTrue(r["limit_exists"]); self.assertEqual(str(r["limit"]), "1")

    def test_algebraic_endpoints(self):
        (lo, hi), = radial_limit("x^2*y/r^3")["components"]      # cos^2 sin: extremes +-2/(3 sqrt 3)
        self.assertAlmostEqual(float(hi), 2 / 27 ** 0.5, places=12)
        self.assertAlmostEqual(float(lo), -2 / 27 ** 0.5, places=12)

    def test_other_point(self):
        self.assertEqual(radial_limit("(x-1)*y/r", point=(1, 0))["set"], "{0}")

    def test_other_square_roots_rejected(self):
        with self.assertRaises(ValueError):
            radial_limit("sqrt(x^2+2*y^2)")

    def test_verify_routes_distance(self):
        r = limit_at_origin("a*b/sqrt(a^2+b^2)")
        self.assertTrue(r["proven"]); self.assertEqual(r["set"], "{0}")


if __name__ == "__main__":
    unittest.main()
