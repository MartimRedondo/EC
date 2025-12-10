import unittest
from random import randint
from polynomials_sage import PolynomialRingDilithium

class TestPolynomialDilithium(unittest.TestCase):
    def setUp(self):
        self.R = PolynomialRingDilithium()

    def test_creation_and_basic_ops(self):
        a = self.R([1] * 256)
        b = self.R([2] * 256)
        self.assertTrue(all(x == 3 for x in (a + b).coeffs))
        self.assertTrue(all(x == 1 for x in (b - a).coeffs))
        self.assertTrue(all(x == 3 for x in (a * 3).coeffs))

    def test_ntt_roundtrip(self):
        a = self.R([i for i in range(256)])
        a_ntt = a.to_ntt()
        a_recovered = a_ntt.from_ntt()
        diff = [(x - y) % self.R.q for x, y in zip(a.coeffs, a_recovered.coeffs)]
        self.assertEqual(sum(diff), 0)

    def test_bit_pack_unpack_t0(self):
        a = self.R([i % 4096 for i in range(256)])
        packed = a.bit_pack_t0()
        unpacked = self.R.bit_unpack_t0(packed)
        self.assertEqual(a.coeffs, unpacked.coeffs)

    def test_sampling_methods(self):
        seed = b"A" * 32
        ball = self.R.sample_in_ball(seed, 39)
        self.assertEqual(len(ball.coeffs), 256)
        self.assertEqual(sum(abs(c) for c in ball.coeffs), 39)

        bounded = self.R.rejection_bounded_poly(seed, 0, eta=2)
        self.assertTrue(all(-2 <= c <= 2 for c in bounded.coeffs))

        poly_ntt = self.R.rejection_sample_ntt_poly(seed, 0, 0)
        self.assertEqual(len(poly_ntt.coeffs), 256)
        self.assertTrue(all(0 <= c < self.R.q for c in poly_ntt.coeffs))

    def test_decompose_and_recompose(self):
        poly = self.R([i for i in range(256)])
        high, low = poly.decompose(alpha=2**10)
        recomposed = high * (2**10) + low
        self.assertEqual(poly.coeffs, recomposed.coeffs)

    def test_hint_operations(self):
        a = self.R([i for i in range(256)])
        b = self.R([(i + 1) for i in range(256)])
        alpha = 2**10
        hint = a.make_hint(b, alpha)
        used = b.use_hint(hint, alpha)
        self.assertEqual(len(used.coeffs), 256)

class TestPolynomialRingDilithium(unittest.TestCase):
    def setUp(self):
        self.R = PolynomialRingDilithium()

    def test_gen(self):
        self.assertEqual(self.R.gen(), self.R([0, 1]))

    def test_non_list_error(self):
        with self.assertRaises(TypeError):
            self.R("not a list")

    def test_long_list_error(self):
        with self.assertRaises(ValueError):
            self.R([0] * (self.R.n + 1))

    def test_string_format(self):
        s = str(self.R)
        self.assertIn("Univariate Polynomial Ring", s)
        self.assertIn("modulus", s)


class TestPolynomialDilithiumElement(unittest.TestCase):
    def setUp(self):
        self.R = PolynomialRingDilithium()

    def test_getitem(self):
        x = self.R.gen()
        self.assertEqual(x[0], 0)
        self.assertEqual(x[1], 1)

    def test_is_zero(self):
        self.assertTrue(self.R(0).is_zero())
        self.assertFalse(self.R(1).is_zero())

    def test_is_constant(self):
        self.assertTrue(self.R(0).is_constant())
        self.assertTrue(self.R(1).is_constant())
        self.assertFalse(self.R.gen().is_constant())

    def test_reduce_coefficients(self):
        for _ in range(100):
            coeffs = [randint(-2 * self.R.q, 3 * self.R.q) for _ in range(self.R.n)]
            f = self.R(coeffs).reduce_coefficients()
            self.assertTrue(all(0 <= c < self.R.q for c in f.coeffs))

    def test_equality_and_repr(self):
        f1 = self.R([1] * self.R.n)
        f2 = self.R([1] * self.R.n)
        self.assertEqual(f1, f2)
        self.assertTrue(str(f1).startswith("1 +"))

    def test_add_sub_mul(self):
        zero = self.R(0)
        one = self.R(1)

        for _ in range(10):
            f1 = self.R.random_element()
            f2 = self.R.random_element()
            f3 = self.R.random_element()

            self.assertEqual(f1 + zero, f1)
            self.assertEqual(f1 * one, f1)
            self.assertEqual(f1 - zero, f1)
            self.assertEqual(f3 - f3, zero)
            self.assertEqual(f1 + f2, f2 + f1)
            self.assertEqual(f1 - f2, -(f2 - f1))
            self.assertEqual(f1 * f2, f2 * f1)
            self.assertEqual(f1 * (f2 * f3), (f1 * f2) * f3)

    def test_pow(self):
        one = self.R(1)
        for _ in range(10):
            f1 = self.R.random_element()
            self.assertEqual(f1 ** 0, one)
            self.assertEqual(f1 ** 1, f1)
            self.assertEqual(f1 * f1, f1 ** 2)
            self.assertEqual(f1 * f1 * f1, f1 ** 3)
            with self.assertRaises(ValueError):
                _ = f1 ** -1

    def test_operator_errors(self):
        f1 = self.R.random_element()
        with self.assertRaises(NotImplementedError):
            _ = f1 + "a"
        with self.assertRaises(NotImplementedError):
            _ = f1 - "a"
        with self.assertRaises(TypeError):
            _ = f1 * "a"
        with self.assertRaises(TypeError):
            _ = f1 ** "a"

    def test_print(self):
        self.assertEqual(str(self.R(0)), "0")
        self.assertEqual(str(self.R(1)), "1")
        self.assertEqual(str(self.R.gen()), "x")
        sample = self.R([1, 2, 3, 4, 1])
        self.assertIn("x^4", str(sample))

if __name__ == "__main__":
    unittest.main()
