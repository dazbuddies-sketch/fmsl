import unittest
from fmsl.optimize import global_min, global_max


class TestGlobalExtrema(unittest.TestCase):
    def test_not_attained(self):
        for f in ("x^2+(x*y-1)^2", "(1-x*y)^2+x^2*y^4", "1/(1+x^2+y^2)"):
            r = global_min(f)
            self.assertEqual(r["value"], "0"); self.assertFalse(r["attained"]); self.assertTrue(r["proven"])

    def test_attained(self):
        r = global_min("x^4+y^4-2*x*y")
        self.assertEqual(r["value"], "-1/2"); self.assertTrue(r["attained"]); self.assertEqual(len(r["points"]), 2)
        r = global_min("x^2*y^4+x^4*y^2-3*x^2*y^2+1")          # Motzkin polynomial
        self.assertEqual(r["value"], "0"); self.assertTrue(r["attained"])
        r = global_max("x/(1+x^2+y^2)")
        self.assertEqual(r["value"], "1/2"); self.assertEqual(r["points"], [("1", "0")])

    def test_unbounded_and_poles(self):
        self.assertEqual(global_min("x^3+y^2")["value"], "-oo")
        with self.assertRaises(ValueError):
            global_min("1/(x^2+y^2-1)")


if __name__ == "__main__":
    unittest.main()
