"""fmsl.reward: problems with proven answers, and partial-credit scoring."""
import unittest

from fmsl.reward import problems, score


class TestReward(unittest.TestCase):
    def test_problems_are_proven_and_varied(self):
        ps = problems(10, seed=3)
        self.assertEqual(len(ps), 10)
        self.assertEqual({p["family"] for p in ps}, {"rational", "analytic", "distance", "infinity", "one_side"})
        for p in ps:
            self.assertIn("ANSWER:", p["prompt"])
            self.assertEqual(score(p, "ANSWER: SET " + p["answer"])["reward"], 1.0)

    def test_partial_credit_orders_answers(self):
        p = {"f": "x^4/(x^2+y^2)^2", "variables": ["x", "y"]}            # truth [0, 1]
        right = score(p, "ANSWER: SET [0, 1]")["reward"]
        close = score(p, "ANSWER: SET [0, 1.2]")["reward"]
        vague = score(p, "ANSWER: SET [-oo, +oo]")["reward"]
        wrong = score(p, "ANSWER: LIMIT 0")["reward"]
        self.assertTrue(right == 1.0 > close > vague > wrong >= 0)

    def test_unparsed(self):
        self.assertEqual(score({"f": "x*y/(x^2+y^2)", "variables": ["x", "y"]}, "no idea")["reward"], 0.0)


if __name__ == "__main__":
    unittest.main()
