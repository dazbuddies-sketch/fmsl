"""sin, cos, exp, log, ... in two and three variables (fmsl.analytic, through fmsl.verify)."""
import unittest

from fmsl.analytic import taylor_quotient, NotAnalytic
from fmsl.verify import limit_at_origin, limit_at, check_claim


class TestAnalytic(unittest.TestCase):
    def test_two_variables(self):
        for e, s in [("sin(x*y)/(x^2+y^2)", "[-1/2, 1/2]"), ("(1-cos(x*y))/(x^2+y^2)^2", "[0, 1/8]"),
                     ("(exp(x)-1-x)/(x^2+y^2)", "[0, 1/2]"), ("log(1+x^2*y)/(x^4+y^2)", "[-1/2, 1/2]"),
                     ("sin(x^2+y^2)/(x^2+y^2)", "{1}"), ("(sqrt(1+x)-1-x/2)/(x^2+y^2)", "[-1/8, 0]"),
                     ("tan(x*y)/(x^2+y^2)", "[-1/2, 1/2]"), ("atan(x*y)/(x^2+y^4)", "[-oo, +oo]"),
                     ("exp(x*y)/(x^2+y^2)", "{+oo}"), ("cos(x)+sin(y)", "{1}")]:
            with self.subTest(e=e):
                r = limit_at_origin(e); self.assertTrue(r["proven"]); self.assertEqual(r["set"], s)

    def test_with_distance(self):
        self.assertEqual(limit_at_origin("sin(x*y)/sqrt(x^2+y^2)")["set"], "{0}")
        self.assertEqual(limit_at_origin("(x-sin(x))/sqrt(x^2+y^2)^3")["set"], "[-1/6, 1/6]")
        self.assertEqual(limit_at_origin("(exp(x*y*z)-1)/(x^2+y^2+z^2)")["set"], "{0}")

    def test_other_point(self):
        self.assertEqual(limit_at("sin((x-1)*y)/((x-1)^2+y^2)", [1, 0])["set"], "[-1/2, 1/2]")
        self.assertEqual(limit_at_origin("x/sqrt((x-1)^2+y^2)")["set"], "{0}")

    def test_common_factor_along_a_line(self):
        for e, st in [("sin(x+y)/(x+y)", "{1}"), ("(1-cos(x+y))/(x+y)^2", "{1/2}"), ("log(1+x*y)/(x*y)", "{1}"),
                      ("tan(2*x-y)/(2*x-y)", "{1}"), ("(x+y)*sin(x*y)/((x+y)*(x^2+y^2))", "[-1/2, 1/2]"),
                      ("sin(x*(x-y))/(x-y)", "{0}"), ("exp(x)*sin(x+y)/(x+y)", "{1}")]:
            with self.subTest(e=e):
                r = limit_at_origin(e); self.assertTrue(r["proven"]); self.assertEqual(r["set"], st)

    def test_blow_up_along_curves(self):
        # the polynomial denominator vanishes along curves the numerator avoids (the curve method)
        for e, st in [("sin(x)/(x-y)", "[-oo, +oo]"), ("sin(x-y)/(x-y)^2", "{-oo} u {+oo}"),
                      ("(1-cos(x))/(x^2-y^2)", "[-oo, 0] u [1/2, +oo]"), ("sin(x)^2/(x^2-y^2)", "[-oo, 0] u [1, +oo]"),
                      ("(x^2+sin(y)^2)/(x*y)", "[-oo, -2] u [2, +oo]"), ("exp(x)/(x-y)", "{-oo} u {+oo}")]:
            with self.subTest(e=e):
                r = limit_at_origin(e); self.assertTrue(r["proven"]); self.assertEqual(r["set"], st)

    def test_curve_zero_not_claimed(self):
        # a common zero curve FMSL cannot cancel, and a denominator that is not a polynomial
        for e in ["(sin(x)-sin(y))/(x-y)", "x/sin(x-y)"]:
            with self.subTest(e=e):
                r = limit_at_origin(e); self.assertFalse(r["proven"]); self.assertIsNone(r["set"])

    def test_taylor_is_exact(self):
        A, B = taylor_quotient("sin(x*y)", ["x", "y"], 7)
        self.assertEqual(A.c, {(1, 1): 1, (3, 3): -__import__("fractions").Fraction(1, 6)})

    def test_refusals(self):
        for e in ["sin(1+x)*y/(x^2+y^2)", "log(2+x)*y/(x^2+y^2)"]:
            with self.subTest(e=e), self.assertRaises(NotAnalytic):
                limit_at_origin(e)

    def test_check_claim(self):
        self.assertTrue(check_claim("sin(x*y)/(x^2+y^2)", "[-1/2, 1/2]")["correct"])
        self.assertFalse(check_claim("sin(x*y)/(x^2+y^2)", "0")["correct"])


if __name__ == "__main__":
    unittest.main()
