"""Three-variable fluid states and the exact two-variable range.  python -m unittest discover tests"""
import unittest
from fractions import Fraction as Fr

from fmsl.space import profile3
from fmsl.plane import image_2d, bgcd, resultant_d
from fmsl.cluster import _union_line, _describe_line


def st(f): return profile3(f)["signed_text"]


class Space(unittest.TestCase):
    def test_limits(self):
        self.assertEqual(profile3("x*y*z/(x^2+y^2+z^2)")["verdict"], "LIMIT = 0 (proven)")
        self.assertTrue(profile3("1/(x^2+y^2+z^2)")["verdict"].startswith("DIVERGES"))

    def test_ranges(self):
        self.assertEqual(st("x*y/(x^2+y^2+z^2)"), "[-1/2, 1/2]")
        self.assertEqual(st("x^2*y/(x^4+y^2+z^2)"), "[-1/2, 1/2]")          # parabola trap in 3-D
        self.assertEqual(st("x^4/((x-y)^2+x^4+z^2)"), "[0, 1]")             # hidden curve, two levels
        self.assertEqual(st("(x*y+y*z+z*x)/(x^2+y^2+z^2)"), "[-1/2, 1]")    # eigenvalues -1/2, 1
        self.assertEqual(st("x^2*y^2*z^2/(x^6+y^6+z^6)"), "[0, 1/3]")       # decided at a vertex
        self.assertEqual(st("(y^2-x*z)/(x^2+z^2+y^4)"), "[-1/2, +oo]")

    def test_irrational_endpoints(self):
        r = profile3("(x*y-2*z^2+x*z)/(x^2+2*y^2+3*z^2)")      # generalized eigenvalues
        lo, hi = r["signed_components"][0]
        self.assertAlmostEqual(float(lo), -0.79687299, places=7)
        self.assertAlmostEqual(float(hi), 0.39497261, places=7)


class Plane(unittest.TestCase):
    def test_curve_zoom(self):
        # the denominator vanishes to high order along the curve yz = x^2 (issue #1)
        self.assertEqual(st("x^4/((y*z-x^2)^2+x^6)"), "[0, +oo]")
        r = profile3("x^2*y*z/((y*z-x^2)^2+x^6)")            # t/(t-1)^2 has minimum -1/4
        self.assertEqual(r["signed_text"], "[-1/4, +oo]")
        self.assertTrue(r["complete"])
        # a common monomial factor must not change the answer
        self.assertEqual(st("x^2*y*z^2/(z*(y*z-x^2)^2+z*x^6)"), "[-1/4, +oo]")

    def test_double_zero_curve_is_not_a_pole_surface(self):
        # (y^2 z - x)^2 + x^4 > 0 off the axes: a squared form touches 0 without crossing it.
        # 0.6.0 wrongly reported [-oo, +oo] here (found by the cube cross-check).
        self.assertEqual(profile3("2*x^3*y^2*z^2/((y^2*z-x)^2+x^4)")["verdict"], "LIMIT = 0 (proven)")
        r = profile3("x^3*y/((y*z-x)^2+x^4)")               # equals 1/z along x = yz
        self.assertEqual(r["signed_text"], "[-oo, +oo]")
        self.assertTrue(r["complete"])

    def test_irrational_zoom(self):
        # zooms at y = +-sqrt(2) x, z = +-sqrt(3) x: arithmetic in Q(sqrt 2, sqrt 3) (issue #2)
        r = profile3("x^4/((y^2-2*x^2)^2+(z^2-3*x^2)^2+x^6)")
        self.assertEqual(r["signed_text"], "[0, +oo]")
        self.assertTrue(r["complete"])
        # x < 0: the denominator vanishes on a real surface; f >= 1/|x| where it is negative
        r = profile3("-x^4/((y^2-3*x^2)^2+x^5+z^2)")
        self.assertEqual(r["signed_text"], "[-oo, 0] u {+oo}")
        self.assertTrue(r["complete"])

    def test_resultant_and_gcd(self):
        P = {(1, 0): Fr(1), (0, 1): Fr(-1)}                 # c - d
        Q = {(2, 0): Fr(1), (0, 2): Fr(1), (0, 0): Fr(-2)}  # c^2 + d^2 - 2
        self.assertEqual(resultant_d(P, Q), [Fr(-2), Fr(0), Fr(2)])     # 2c^2 - 2
        G = bgcd({(2, 0): Fr(1), (0, 2): Fr(-1)}, {(1, 0): Fr(1), (0, 1): Fr(1)})
        self.assertEqual(len(G), 2)                         # c + d

    def test_range_over_plane_and_box(self):
        cd = {(1, 1): Fr(1)}; den = {(0, 0): Fr(1), (4, 0): Fr(1), (0, 2): Fr(1)}
        iv, _ = image_2d(cd, den)                           # cd/(1+c^4+d^2): [-sqrt2/4, sqrt2/4]
        lo, hi = _union_line(iv)[0]
        self.assertAlmostEqual(float(hi), 2 ** 0.5 / 4, places=9)
        iv, ex = image_2d({(0, 0): Fr(1)}, {(2, 2): Fr(3)})  # 1/(3 c^2 d^2): [0, +oo]
        self.assertEqual(_describe_line(_union_line(iv)), "[0, +oo]")


if __name__ == "__main__":
    unittest.main()
