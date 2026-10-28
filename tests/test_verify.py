import math, unittest
from fmsl.verify import limit_at_origin, check_claim, parse_claim


class TestLimitAtOrigin(unittest.TestCase):
    def test_one_variable(self):
        self.assertEqual(limit_at_origin("sin(x)/x")["limit"], "1")
        self.assertEqual(limit_at_origin("(1-cos(t))/t^2")["limit"], "1/2")
        r = limit_at_origin("1/x")
        self.assertFalse(r["limit_exists"]); self.assertEqual(r["set"], "{-oo} u {+oo}")
        self.assertEqual(limit_at_origin("x*sin(1/x)")["limit"], "0")
        self.assertEqual(limit_at_origin("abs(x)/x")["set"], "{-1} u {1}")

    def test_two_and_three_variables(self):
        r = limit_at_origin("x^3*y/(x^6+y^2)")
        self.assertEqual(r["set"], "[-1/2, 1/2]"); self.assertTrue(r["proven"]); self.assertTrue(r["witnesses"])
        self.assertEqual(limit_at_origin("(a^2+b^2)/(a*b)")["set"], "[-oo, -2] u [2, +oo]")
        self.assertEqual(limit_at_origin("x*y*z/(x^2+y^2+z^2)")["limit"], "0")
        self.assertEqual(limit_at_origin("(x*y+y*z+z*x)/(x^2+y^2+z^2)")["set"], "[-1/2, 1]")


class TestCheckClaim(unittest.TestCase):
    def test_parse(self):
        self.assertEqual(parse_claim("LIMIT 1/2"), [(0.5, 0.5)])
        self.assertEqual(parse_claim("[-oo, 1-2sqrt(3)] u {0}")[0][0], -math.inf)
        self.assertEqual(parse_claim("{-oo, +oo}"), [(-math.inf, -math.inf), (math.inf, math.inf)])

    def test_verdicts(self):
        f = "-2*x^2/(x^2+2*x*y-2*y^2+x^3*y)"
        r = check_claim(f, "[-oo, +oo]")
        self.assertFalse(r["correct"]); self.assertEqual(r["extra"], "(-4/3, 0)"); self.assertIsNone(r["missing"])
        self.assertTrue(check_claim(f, "[-oo, -4/3] u [0, +oo]")["correct"])
        r = check_claim("x^3*y/(x^6+y^2)", "LIMIT 0")
        self.assertFalse(r["correct"]); self.assertEqual(r["missing"], "[-1/2, 0) u (0, 1/2]")
        self.assertTrue(check_claim("sin(x)/x", "1")["correct"])
        self.assertFalse(check_claim("1/x", "+oo")["correct"])


if __name__ == "__main__":
    unittest.main()


class TestOtherPoints(unittest.TestCase):
    def test_point_and_side(self):
        from fmsl.verify import limit_at
        self.assertEqual(limit_at("(x^2-4)/(x-2)", 2)["limit"], "4")
        self.assertEqual(limit_at("1/(x-2)", 2, side="-")["one_sided_limit"], "-oo")
        self.assertEqual(limit_at("abs(x-2)/(x-2)", 2, side="+")["one_sided_limit"], "1")
        self.assertEqual(limit_at("(x*y-2*y)/((x-2)^2+y^2)", (2, 0))["set"], "[-1/2, 1/2]")
