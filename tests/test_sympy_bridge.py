import unittest

try:
    import sympy as sp
except ImportError:
    sp = None


@unittest.skipIf(sp is None, "SymPy not installed")
class TestSympyBridge(unittest.TestCase):
    def setUp(self):
        self.x, self.y, self.z = sp.symbols("x y z", real=True)

    def test_limit_set(self):
        from fmsl.sympy_bridge import limit_set
        x, y = self.x, self.y
        self.assertEqual(limit_set(x**3*y/(x**6 + y**2), (x, y), (0, 0)), sp.Interval(-sp.Rational(1, 2), sp.Rational(1, 2)))
        s = limit_set((x**2 + y**2)/(x*y), (x, y))
        self.assertTrue(sp.oo in s and -sp.oo in s and 2 in s and 0 not in s)
        # an algebraic endpoint (DGNP 2022, example 4): [-sqrt(2)/4, sqrt(2)/4]
        f = (4*x**2*y**2 - 4*x*y**3 + y**4 - 2*x*y**2 + y**3)/(8*x**2*y**2 - 8*x*y**3 + 3*y**4 + 8*x**2 - 8*x*y + 2*y**2)
        s = limit_set(f, (x, y))
        self.assertAlmostEqual(float(s.sup), 2**0.5/4, places=12)
        self.assertEqual(sp.nsimplify(s.sup), sp.sqrt(2)/4)

    def test_limit_and_point(self):
        from fmsl.sympy_bridge import multivariate_limit, NoLimit
        x, y, z = self.x, self.y, self.z
        self.assertEqual(multivariate_limit(x**2*y**2/(x**2 + y**2), (x, y)), 0)
        # at the point (1, 2): (x-1)(y-2)^2/((x-1)^2+(y-2)^2) -> 0
        self.assertEqual(multivariate_limit((x - 1)*(y - 2)**2/((x - 1)**2 + (y - 2)**2), (x, y), (1, 2)), 0)
        with self.assertRaises(NoLimit) as cm:
            multivariate_limit((x*y + y*z + z*x)/(x**2 + y**2 + z**2), (x, y, z))
        self.assertEqual(cm.exception.values, sp.Interval(-sp.Rational(1, 2), 1))


if __name__ == "__main__":
    unittest.main()
