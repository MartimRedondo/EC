import unittest
from modules_sage import ModuleDilithium, MatrixDilithium
from polynomials_sage import PolynomialRingDilithium
import random

class TestMatrixDilithium(unittest.TestCase):

    def setUp(self):
        self.R = PolynomialRingDilithium()
        self.A = MatrixDilithium.from_seed(b"A"*32, 2, 2, self.R)
        self.B = MatrixDilithium.from_seed(b"B"*32, 2, 2, self.R)

    def test_addition(self):
        C = self.A + self.B
        for i in range(2):
            for j in range(2):
                self.assertEqual(C[i, j], self.A[i, j] + self.B[i, j])

    def test_subtraction(self):
        C = self.A - self.B
        for i in range(2):
            for j in range(2):
                self.assertEqual(C[i, j], self.A[i, j] - self.B[i, j])

    def test_multiplication(self):
        A = MatrixDilithium.from_seed(b"A"*32, 2, 3, self.R)
        B = MatrixDilithium.from_seed(b"B"*32, 3, 2, self.R)
        C = A * B
        self.assertEqual(C.shape(), (2, 2))

    def test_transpose(self):
        T = self.A.transpose()
        for i in range(2):
            for j in range(2):
                self.assertEqual(self.A[i, j], T[j, i])

    def test_map(self):
        mapped = self.A.map(lambda x: x * 2)
        for i in range(2):
            for j in range(2):
                self.assertEqual(mapped[i, j], self.A[i, j] * 2)

    def test_ntt_conversion(self):
        self.ring = PolynomialRingDilithium()
        self.A = MatrixDilithium([
            [self.ring([1] * 256), self.ring([2] * 256)],
            [self.ring([3] * 256), self.ring([4] * 256)]
        ], self.ring)

        A_ntt = self.A.to_ntt()
        A_back = A_ntt.from_ntt()

        for i in range(2):
            for j in range(2):
                self.assertEqual(self.A[i, j].coeffs, A_back[i, j].coeffs)

    def test_sample_vector_ntt(self):
        v = self.A.sample_vector_ntt(b"A"*32, 3)
        self.assertEqual(v.shape(), (3, 1))


class TestModuleDilithium(unittest.TestCase):

    def setUp(self):
        self.module = ModuleDilithium()
        self.R = self.module.ring

    def test_bit_unpack_t0(self):
        packed = bytes([random.randint(0, 255) for _ in range(416 * 2)])
        M = self.module.bit_unpack_t0(packed, 2, 1)
        self.assertEqual(M.shape(), (2, 1))

    def test_bit_unpack_t1(self):
        packed = bytes([random.randint(0, 255) for _ in range(320 * 2)])
        M = self.module.bit_unpack_t1(packed, 2, 1)
        self.assertEqual(M.shape(), (2, 1))

    def test_bit_unpack_s(self):
        packed = bytes([random.randint(0, 255) for _ in range(96 * 2)])
        M = self.module.bit_unpack_s(packed, 2, 1, eta=2)
        self.assertEqual(M.shape(), (2, 1))

    def test_bit_unpack_w(self):
        packed = bytes([random.randint(0, 255) for _ in range(192 * 2)])
        M = self.module.bit_unpack_w(packed, 2, 1, gamma_2=95232)
        self.assertEqual(M.shape(), (2, 1))

    def test_bit_unpack_z(self):
        packed = bytes([random.randint(0, 255) for _ in range(576 * 2)])
        M = self.module.bit_unpack_z(packed, 2, 1, gamma_1=(1 << 17))
        self.assertEqual(M.shape(), (2, 1))


if __name__ == "__main__":
    unittest.main()
