"""Run with:  python -m unittest discover tests"""
import random
import unittest
from fractions import Fraction as Fr

from fmsl import V, x, Exp, Sin, Cos, Essential, compile_text, I, proven_series_error
from fmsl.compiler import Parser, build, series_data, series_eval, ref_eval, rel


def lead(v):
    r = v.val()
    if r[0] == "zero":
        return "C_e" if r[1] else "unproven"
    return (r[1], v.coeff(r[1]))


class ClassicLimits(unittest.TestCase):
    def test_removable(self):
        self.assertEqual(lead(Sin(x) / x), (0, 1))
        self.assertEqual(lead((1 - Cos(x)) / x**2), (0, Fr(1, 2)))
        self.assertEqual(lead((Exp(x) - 1 - x) / x**2), (0, Fr(1, 2)))
        self.assertEqual(lead((x - Sin(x)) / x**3), (0, Fr(1, 6)))

    def test_pole(self):
        self.assertEqual(lead((Sin(x) - x) / x**4), (-1, Fr(-1, 6)))
        self.assertEqual(lead(1 / x), (-1, 1))

    def test_zero_over_zero_states_differ(self):
        self.assertEqual(lead(x / x), (0, 1))
        self.assertEqual(lead(x**2 / x), (1, 1))

    def test_rearrangement_is_consistent(self):
        f, g = x, 1 / x
        self.assertEqual(lead((f + g) - g), lead(f + (g - g)))


class CertifiedZero(unittest.TestCase):
    def test_identities(self):
        for h in [Sin(x)**2 + Cos(x)**2 - 1,
                  Sin(2 * x) - 2 * Sin(x) * Cos(x),
                  Exp(2 * x) - Exp(x)**2,
                  Cos(3 * x) - 4 * Cos(x)**3 + 3 * Cos(x),
                  Exp(x**2) * Exp(-x**2) - 1,
                  Sin(x**3)**2 + Cos(x**3)**2 - 1,
                  Cos(x + x**2) - Cos(x) * Cos(x**2) + Sin(x) * Sin(x**2)]:
            self.assertEqual(lead(h), "C_e")

    def test_z_minus_z(self):
        z = Exp(x) / x
        self.assertEqual(lead(z - z), "C_e")

    def test_near_misses_are_not_zero(self):
        self.assertEqual(lead(Sin(x)**2 + Cos(x)**2 - 1 + x**9 / 1000), (9, Fr(1, 1000)))
        self.assertEqual(lead(Exp(x**2) * Exp(-x**2) - 1 + x**20 / 7), (20, Fr(1, 7)))
        self.assertEqual(lead(Exp(x**61) - 1), (61, 1))

    def test_high_order_vanishing(self):
        self.assertEqual(lead((Sin(x) - x)**10)[0], 30)
        self.assertEqual(lead((Exp(x**2) - 1 - x**2)**5)[0], 20)
        self.assertEqual(lead(Exp(x**40) - 1 - x**40)[0], 80)

    def test_random_identities_never_wrong(self):
        rng = random.Random(1)
        for _ in range(40):
            a = Fr(rng.choice([1, 2, -1, 3]), rng.choice([1, 2, 3]))
            h = Sin(2 * a * x) - 2 * Sin(a * x) * Cos(a * x)
            if rng.random() < 0.5:
                k = rng.randint(1, 12)
                self.assertEqual(lead(h + x**k), (k, 1))
            else:
                self.assertEqual(lead(h), "C_e")


class Boundaries(unittest.TestCase):
    def test_uncertified_is_reported_not_guessed(self):
        self.assertEqual(lead(Exp(Sin(x)) * Exp(-Sin(x)) - 1), "unproven")

    def test_essential(self):
        with self.assertRaises(Essential):
            Exp(1 / x**2)
        with self.assertRaises(Essential):
            Sin(1 / x)

    def test_transcendental_constant_refused(self):
        with self.assertRaises(ValueError):
            Exp(1 + x)


class Language(unittest.TestCase):
    def test_program(self):
        out = "\n".join(compile_text(
            "fluid B(t) = (1 - cos(t))/t^2\n"
            "at t -> 0: B\n"
            "check sin(t)^2 + cos(t)^2 - 1\n"
            "check exp(sin(t))*exp(-sin(t)) - 1\n"
            "if B == 0.5\n"
            "emit python B\n"))
        self.assertIn("removable: value 1/2 [proven]", out)
        self.assertIn("proven_zero", out)
        self.assertIn("unproven", out)
        self.assertIn("TYPE ERROR", out)
        self.assertIn("def B(t):", out)

    def test_emitted_python_is_accurate(self):
        out = "\n".join(compile_text("fluid B(t) = (1 - cos(t))/t^2\nemit python B\n"))
        code = out.split("emit python B:\n", 1)[1].split("        accuracy")[0]
        code = "\n".join(line[8:] for line in code.splitlines())
        ns = {}
        exec(code, ns)
        for t in [0.0, 1e-12, 1e-8, -1e-5, 0.3, 2.0]:
            exact = 0.5 if t == 0 else None
            if exact is not None:
                self.assertEqual(ns["B"](t), exact)
        self.assertAlmostEqual(ns["B"](1e-8), 0.5, places=15)
        self.assertAlmostEqual(ns["B"](-1e-5), 0.5 - 1e-10 / 24, places=15)


class ComplexValues(unittest.TestCase):
    def test_euler_is_proven(self):
        self.assertEqual(lead(Exp(I * x) - Cos(x) - I * Sin(x)), "C_e")
        self.assertEqual(lead(Exp(I * x) * Exp(-I * x) - 1), "C_e")

    def test_quantum_limits(self):
        self.assertEqual(lead((Exp(-3 * I * x) - 1) / (-I * x)), (0, 3))
        self.assertEqual(lead((Exp(I * x) - 1 - I * x) / (I * x)**2), (0, Fr(1, 2)))
        self.assertEqual(lead((Exp(-3 * I * x) - 1 + 3 * I * x) / x**2), (0, Fr(-9, 2)))

    def test_near_miss_complex(self):
        self.assertEqual(lead(Exp(I * x) - Cos(x) - I * Sin(x) + x**9), (9, 1))

    def test_infinity_flow_direction(self):
        text = (I / x).describe()
        self.assertIn("90.0 deg from the right", text)
        self.assertIn("-90.0 deg from the left", text)

    def test_spinning_phase_is_essential(self):
        with self.assertRaises(Essential):
            Exp(I / x)


class SiletzProofs(unittest.TestCase):
    def test_bound_is_never_below_actual_error(self):
        h0 = Fr(1, 2)
        pts = [h0 * Fr(k, 97) for k in range(-96, 97) if k]
        for src, var, cx in [("(exp(i*z) - 1 - i*z)/(i*z)^2", "z", True),
                             ("(exp(-3*g) - 1 + 3*g)/g^2", "g", False),
                             ("exp(x)/x^2 + 1/x", "x", False)]:
            v = build(Parser(src).parse(), var)
            lo, co, _ = series_data(v)
            for K in (4, 8, 14):
                pb, _ = proven_series_error(v, lo, co, K, h0, cx, abs(lo))
                if pb is None:
                    continue
                meas = max(rel(series_eval(lo, co, K, float(p)), ref_eval(lo, co, float(p))) for p in pts)
                self.assertLessEqual(meas, float(pb), f"{src} K={K}")

    def test_emit_carries_proof(self):
        out = "\n".join(compile_text("fluid B(t) = (1 - cos(t))/t^2\nemit python B\n"))
        self.assertIn("Siletz proof: |rel. error| <=", out)


class Context(unittest.TestCase):
    def test_context_parameter(self):
        out = "\n".join(compile_text("context T = 3\nfluid p(w) = (exp(-T*i*w) - 1)/(-i*w)\nat w -> 0: p\n"))
        self.assertIn("removable: value 3 [proven]", out)

    def test_context_must_be_exact(self):
        out = "\n".join(compile_text("context T = sin(1)\n"))
        self.assertIn("ERROR", out)


if __name__ == "__main__":
    unittest.main()
