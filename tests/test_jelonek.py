import unittest
from fmsl.jelonek import jelonek_set


class TestJelonek(unittest.TestCase):
    def test_known(self):
        self.assertEqual(jelonek_set("x", "x*y")["text"], "the line a = 0")
        self.assertTrue(jelonek_set("x+y^2", "y")["empty"])              # an automorphism
        self.assertTrue(jelonek_set("x", "y")["empty"])
        r = jelonek_set("x*y", "x*y^2-y")
        self.assertEqual(sorted((L["fixed"], str(L["value"])) for L in r["lines"]), [("a", "1"), ("b", "0")])
        self.assertEqual(jelonek_set("x^2*y^2+x*y", "y")["text"], "the line b = 0, a in [-1/4, +oo)")
        r = jelonek_set("x*y", "x^2*y^2+x")
        self.assertEqual(len(r["curves"]), 1); self.assertEqual(r["curves"][0]["B"][-1], 1)   # b = a^2


if __name__ == "__main__":
    unittest.main()
