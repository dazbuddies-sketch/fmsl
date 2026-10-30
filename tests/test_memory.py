"""fmsl.memory: saved proofs and strategy learning (both off by default)."""
import os, tempfile, unittest

import fmsl.memory as memory
from fmsl.verify import limit_at_origin


class TestMemory(unittest.TestCase):
    def tearDown(self):
        memory.disable()

    def test_off_by_default(self):
        self.assertFalse(memory.enabled())
        self.assertEqual(memory.order("k", 3), [0, 1, 2])

    def test_saved_answer_survives_a_new_session(self):
        path = os.path.join(tempfile.mkdtemp(), "m.json")
        memory.enable(path)
        r1 = limit_at_origin("x*y/(x^2+y^2)")
        memory.disable(); memory.enable(path)
        r2 = limit_at_origin("x*y/(x^2+y^2)")
        self.assertEqual(r1["set"], r2["set"])
        self.assertTrue(r2["explanation"].startswith("from memory"))

    def test_unproven_answers_are_not_saved(self):
        path = os.path.join(tempfile.mkdtemp(), "m.json")
        memory.enable(path)
        limit_at_origin("x/sin(x-y)")
        self.assertEqual(memory.stats()["results"], 0)

    def test_failure_within_a_time_limit_is_remembered(self):
        import time
        path = os.path.join(tempfile.mkdtemp(), "m.json")
        memory.enable(path)
        f = "x*y*z/sqrt(x^2+y^2+z^2)^3"                       # takes well over 2 s to prove
        self.assertFalse(limit_at_origin(f, time_limit=2)["proven"])
        memory.disable(); memory.enable(path)
        t = time.time(); r = limit_at_origin(f, time_limit=2)
        self.assertFalse(r["proven"]); self.assertLess(time.time() - t, 0.5)
        self.assertIn("earlier attempt", r["explanation"])

    def test_bandit_prefers_what_worked(self):
        memory.enable(os.path.join(tempfile.mkdtemp(), "m.json"))
        for _ in range(5):
            memory.reward("kind", 0, 0.0); memory.reward("kind", 1, 0.9); memory.reward("kind", 2, 0.2)
        self.assertEqual(memory.order("kind", 3)[0], 1)
        self.assertEqual(memory.order("other kind", 3), [0, 1, 2])     # nothing learned there yet

    def test_new_version_keeps_strategy_only(self):
        path = os.path.join(tempfile.mkdtemp(), "m.json")
        memory.enable(path); limit_at_origin("x*y/(x^2+y^2)"); memory.reward("kind", 1, 1.0)
        data = memory._STATE["data"]; data["version"] = "0.0.1"; memory._save()
        memory.disable(); memory.enable(path)
        self.assertEqual(memory.stats()["results"], 0)
        self.assertIn("kind", memory.stats()["strategy"])


class TestPortfolio(unittest.TestCase):
    def test_another_start_proves_what_the_first_cannot(self):
        # the first coordinate direction is not proven here; another start is, in well under a second
        f = ("((-1)*x^2*y^0*z^0*sqrt(x^2+y^2+z^2)^1+(-2)*x^0*y^2*z^0*sqrt(x^2+y^2+z^2)^0)/"
             "((1)*x^1*y^0*z^1*sqrt(x^2+y^2+z^2)^0)")
        r = limit_at_origin(f, ["x", "y", "z"])
        self.assertTrue(r["proven"])


if __name__ == "__main__":
    unittest.main()
