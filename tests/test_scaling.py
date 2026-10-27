import unittest
from fractions import Fraction as Fr

from fmsl.multivar import parse_rational
from fmsl.scaling import balance_scales, regimes, scaling_map


class ScalingMap(unittest.TestCase):
    def test_balance_scales(self):
        N, D = parse_rational("x^2*y/(x^4+y^2)")
        self.assertEqual(balance_scales(N, D), [Fr(2)])
        N, D = parse_rational("(x^2+y^3)/(x^2+y^2)")
        self.assertEqual(balance_scales(N, D), [Fr(2, 3), Fr(1)])

    def test_regime_values(self):
        N, D = parse_rational("(x^2-y^2)/(x^2+y^2)")
        beh = [g["behaviour"] for g in regimes(N, D) if g["kind"] == "regime"]
        self.assertEqual(beh, ["f -> -1", "f -> 1"])

    def test_hidden_scale_is_flagged(self):
        N, D = parse_rational("x^4/((x-y)^2+x^4)")
        bal = [g for g in regimes(N, D) if g["kind"] == "balance"]
        self.assertIn("a finer scale decides there", bal[0]["behaviour"])
        r = scaling_map("x^4/((x-y)^2+x^4)")
        ends = dict((str(v), p) for v, p in r["deciding"])
        self.assertTrue(any("exactly" in p for p in ends["1"]))     # the end 1 lives on y = x

    def test_ends_decided_at_balance(self):
        r = scaling_map("x^2*y/(x^4+y^2)")
        for v, paths in r["deciding"]:
            self.assertTrue(all("|x|^2" in p for p in paths))


if __name__ == "__main__":
    unittest.main()
