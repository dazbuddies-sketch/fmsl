"""Limits at infinity and from one side in several variables (fmsl.regions, through fmsl.verify)."""
import unittest

from fmsl.verify import limit_at
from fmsl.regions import RegionError


class TestRegions(unittest.TestCase):
    def test_at_infinity(self):
        for e, s in [("x*y/(x^2+y^2)", "[-1/2, 1/2]"), ("(x^2+y)/(x^2+y^2)", "[0, 1]"),
                     ("(x^2-y^2)/(x^2+y^2+1)", "[-1, 1]"), ("x/(x^2+y^2)", "{0}"),
                     ("x*y/sqrt(x^2+y^2)", "[-oo, +oo]"), ("sqrt(x^2+2*y^2)/(x+2*y^2)", "[-oo, -1] u [0, +oo]")]:
            with self.subTest(e=e):
                r = limit_at(e, "oo"); self.assertTrue(r["proven"]); self.assertEqual(r["set"], s)

    def test_one_sided(self):
        for e, sd, s in [("x*y/(x^2+y^2)", {"x": "+", "y": "+"}, "[0, 1/2]"),
                         ("x*y/(x^2+y^2)", {"x": "+", "y": "-"}, "[-1/2, 0]"),
                         ("y/x", {"x": "+", "y": "+"}, "[0, +oo]"),
                         ("sin(x*y)/(x^2+y^2)", {"x": "+"}, "[-1/2, 1/2]")]:
            with self.subTest(e=e, sd=sd):
                r = limit_at(e, [0, 0], side=sd); self.assertTrue(r["proven"]); self.assertEqual(r["set"], s)

    def test_other_point(self):
        self.assertEqual(limit_at("(x-1)*y/((x-1)^2+y^2)", [1, 0], side={"x": "-", "y": "+"})["set"], "[-1/2, 0]")   # x - 1 < 0 < y

    def test_refusals(self):
        with self.assertRaises(RegionError):
            limit_at("sin(x*y)/(x^2+y^2)", "oo")
        with self.assertRaises(RegionError):
            limit_at("x/sqrt(x^2+y^2)", [0, 0], side={"x": "+"})


if __name__ == "__main__":
    unittest.main()
