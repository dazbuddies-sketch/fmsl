"""Two-variable fluid states.  python -m unittest discover tests"""
import unittest
from fractions import Fraction as Fr

from fmsl.multivar import profile, _ueval


def verdict(f):
    return profile(f)["verdict"][0]


def f_at(N, D, x, y):
    d = sum(v * x ** i * y ** j for (i, j), v in D.items())
    return None if d == 0 else sum(v * x ** i * y ** j for (i, j), v in N.items()) / d


class Classic(unittest.TestCase):
    def test_limits_that_exist(self):
        for f in ["x^3/(x^2+y^2)", "x^2*y^2/(x^2+y^2)", "(x^3+y^3)/(x^2+y^2)", "(x^2-y^2)/(x-y)"]:
            self.assertEqual(verdict(f), "LIMIT = 0 (proven)", f)
        self.assertEqual(verdict("(x^2+y^2)/(x^2+y^2)"), "LIMIT = 1 (proven)")

    def test_no_limit(self):
        for f in ["x*y/(x^2+y^2)", "x/y", "x*y^2/(x^2+y^4)"]:
            self.assertEqual(verdict(f), "NO LIMIT (proven)", f)

    def test_divergence(self):
        self.assertTrue(verdict("1/(x^2+y^2)").startswith("DIVERGES"))


class HiddenCurves(unittest.TestCase):
    def test_parabola_trap(self):
        r = profile("x^2*y/(x^4+y^2)")
        self.assertEqual(r["verdict"][0], "NO LIMIT (proven)")
        lines = [e for e in r["entries"] if e.get("a") == 1 and e["kind"] == "value"]
        self.assertTrue(lines and all(e["value"] == 0 for e in lines))          # every line -> 0
        dep = [e for e in r["entries"] if e["kind"] == "depends"]
        self.assertEqual(dep[0]["a"], 2)                                        # balance at a = 2
        self.assertEqual(_ueval(dep[0]["R"][0], Fr(1)) / _ueval(dep[0]["R"][1], Fr(1)), Fr(1, 2))

    def test_two_levels_deep(self):
        r = profile("x^4/(x^2 + y^2 - 2*x*y + x^4)")
        self.assertEqual(r["verdict"][0], "NO LIMIT (proven)")
        self.assertTrue(any(e["kind"] == "depends" and e.get("depth") == 1 for e in r["entries"]))

    # --- regressions for bugs found during development ---
    def test_denominator_root_with_positive_order_gap(self):
        # (x - y)^2 in the denominator: blows up near y = x although the numerator 'looks' smaller
        self.assertEqual(verdict("x^2*y/(x^2 + y^2 - 2*x*y)"), "NO LIMIT (proven)")

    def test_hidden_zero_curve_of_denominator(self):
        # for x < 0 the denominator vanishes on y = x^3 +- |x|^3.5; f is unbounded beside it
        self.assertEqual(verdict("(y-x^3)^2*x/((y-x^3)^2 + x^7)"), "NO LIMIT (proven)")

    def test_cross_balance_between_numerator_and_denominator(self):
        # y^2/(x(2x-1)): blows up along y = x^(1/3); the switch happens where N ties D (a = 1/2)
        self.assertEqual(verdict("((-1)*x^3*y^2)/((-1)*x^4 + (2)*x^5)"), "NO LIMIT (proven)")

    def test_pole_curve_by_implicit_function_theorem(self):
        # (-2-3x^3)/(2x^4 (x - xy + y^3)): the pole curve x = -y^3/(1-y) has an infinite Puiseux
        # series; a simple root of the denominator's form settles it without zooming forever
        f = "((-2)*x^1*y^0 + (-3)*x^4*y^0)/((2)*x^6*y^0 + (-2)*x^6*y^1 + (2)*x^5*y^3)"
        self.assertTrue(verdict(f).startswith("DIVERGES"))

    def test_linear_root_over_algebraic_field(self):
        f = "((-1)*x^2*y^0 + (1)*x^3*y^1 + (3)*x^0*y^5)/((-2)*x^4*y^6 + (3)*x^4*y^0 + (3)*x^2*y^6)"
        self.assertEqual(verdict(f), "NO LIMIT (proven)")

    def test_deep_path_exponents_are_in_x(self):
        r = profile("x^6/((y-x-x^2)^2+x^6)")                       # three zoom levels
        self.assertIn("y = 1*|x|^1  (x>0) + 1*|x|^2 + c*|x|^3", r["verdict"][1][0])

    def test_every_value_beside_simple_curves(self):
        # x^3/(y - x^2): along y = x^2 + c x^3 the value is 1/c -> every real value
        self.assertIn("every real value", " ".join(profile("x^3/(y-x^2)")["verdict"][1]))
        self.assertTrue(verdict("1/y").startswith("DIVERGES"))     # pole curve, but infinite all around


class DirectEvaluation(unittest.TestCase):
    def test_every_reported_path_matches_direct_evaluation(self):
        for f in ["x^2*y/(x^4+y^2)", "x^4/(x^2 + y^2 - 2*x*y + x^4)", "x*y/(x^2+y^2)", "x^3/(x^2+y^2)"]:
            r = profile(f); N, D = r["N"], r["D"]
            for e in r["entries"]:
                if "point" not in e or e["kind"] not in ("value", "depends"):
                    continue
                for c in ((e["fixed_c"],) if "fixed_c" in e else (Fr(1, 3), Fr(5, 2))):
                    if e["kind"] == "depends":
                        dn = _ueval(e["R"][1], c)
                        if dn == 0: continue
                        pred = _ueval(e["R"][0], c) / dn
                    else:
                        pred = e["value"]
                    x, y = e["point"](c, Fr(1, 10 ** 15))
                    val = f_at(N, D, x, y)
                    if val is None: continue
                    self.assertLess(abs(val - pred), Fr(1, 10 ** 4) * (1 + abs(pred)), f"{f}: {e['path']} c={c}")


if __name__ == "__main__":
    unittest.main()
