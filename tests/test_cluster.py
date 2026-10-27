"""Exact fluid-state arcs (cluster sets).  python -m unittest discover tests"""
import unittest
from fractions import Fraction as Fr

from fmsl.multivar import profile
from fmsl.cluster import image_of_rational_function, POS, NEG, cmp
from fmsl.algebraic import A


def arc(f): return profile(f)["cluster"]


class Arcs(unittest.TestCase):
    def test_points(self):
        self.assertEqual(arc("x^3/(x^2+y^2)")["text"], "{0}")
        self.assertEqual(arc("1/(x^2+y^2)")["text"], "{oo}")

    def test_parabola_trap_is_exactly_half_interval(self):
        c = arc("x^2*y/(x^4+y^2)")
        self.assertEqual(c["text"], "[-1/2, 1/2]")
        self.assertTrue(c["exact"] and c["connected"])

    def test_hidden_two_level_curve(self):
        self.assertEqual(arc("x^4/((x-y)^2+x^4)")["text"], "[0, 1]")

    def test_arc_through_infinity(self):
        # c/(1 - 3c + c^2): poles at (3 +- sqrt5)/2, critical values -1 (c = 1) and -1/5 (c = -1)
        c = arc("x*y/(x^2-3*x*y+y^2)")
        self.assertTrue(c["connected"] and not c["whole"])
        self.assertIn("from -1/5 up through oo to -1", c["text"])

    def test_algebraic_endpoint(self):
        c = arc("x^2/(x^2+2*y^2)+y^2/(3*x^2+y^2)")
        lo = c["components"][0][0]
        self.assertIsInstance(lo, A)
        self.assertTrue(Fr(5797958, 10 ** 7) < lo < Fr(5797959, 10 ** 7))
        self.assertTrue(c["exact"])

    def test_every_value(self):
        self.assertTrue(arc("x/y")["whole"])
        self.assertTrue(arc("x^3/(y-x^2)")["whole"])

    def test_image_of_rational_function(self):
        iv = image_of_rational_function([0, 1], [1, 0, 1])         # c/(1+c^2)
        lo = min((p[0] for p in iv), key=lambda v: float(v)); hi = max((p[1] for p in iv), key=lambda v: float(v))
        self.assertEqual((lo, hi), (Fr(-1, 2), Fr(1, 2)))


if __name__ == "__main__":
    unittest.main()
