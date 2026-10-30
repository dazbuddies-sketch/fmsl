"""fmsl.indeterminacy: clean indeterminacy lines in three variables force every value."""
import unittest

from fmsl.indeterminacy import clean_line
from fmsl.space import parse3
from fmsl.verify import limit_at_origin


class TestIndeterminacy(unittest.TestCase):
    def test_clean_lines_found(self):
        for e in ["(x*y+z)/(x*z-y)", "(x+y*z)/(x-y*z)", "x/(y+z)"]:
            with self.subTest(e=e):
                self.assertIsNotNone(clean_line(*parse3(e)))
                r = limit_at_origin(e, ["x", "y", "z"])
                self.assertTrue(r["proven"]); self.assertEqual(r["set"], "[-oo, +oo]")
                self.assertIn("indeterminacy", r["explanation"])

    def test_not_clean_not_claimed(self):
        # the z-axis is in {P = Q = 0} but dP = 0 there: the test must not fire, and the answer is [0, 1]
        self.assertIsNone(clean_line(*parse3("x^2/(x^2+y^2)+0*z")))
        self.assertEqual(limit_at_origin("x^2/(x^2+y^2+0*z)", ["x", "y", "z"])["set"], "[0, 1]")


if __name__ == "__main__":
    unittest.main()
