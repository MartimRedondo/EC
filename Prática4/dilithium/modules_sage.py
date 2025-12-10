from typing import List
from sage.all import *
import hashlib

from polynomials_sage import PolynomialRingDilithium 

class ModuleDilithium:
    def __init__(self):
        self.ring = PolynomialRingDilithium()
        self.matrix = MatrixDilithium

    def __bit_unpack(self, input_bytes, m, n, alg, packed_len, *args):
        poly_bytes = [
            input_bytes[i : i + packed_len]
            for i in range(0, len(input_bytes), packed_len)
        ]
        matrix = [
            [alg(poly_bytes[n * i + j], *args) for j in range(n)] for i in range(m)
        ]
        return MatrixDilithium(matrix, self.ring)

    def bit_unpack_t0(self, input_bytes, m, n):
        return self.__bit_unpack(input_bytes, m, n, self.ring.bit_unpack_t0, 416)

    def bit_unpack_t1(self, input_bytes, m, n):
        return self.__bit_unpack(input_bytes, m, n, self.ring.bit_unpack_t1, 320)

    def bit_unpack_s(self, input_bytes, m, n, eta):
        if eta == 2:
            packed_len = 96
        elif eta == 4:
            packed_len = 128
        else:
            raise ValueError("eta inválido")
        return self.__bit_unpack(input_bytes, m, n, self.ring.bit_unpack_s, packed_len, eta)

    def bit_unpack_w(self, input_bytes, m, n, gamma_2):
        if gamma_2 == 95232:
            packed_len = 192
        elif gamma_2 == 261888:
            packed_len = 128
        else:
            raise ValueError("gamma_2 inválido")
        return self.__bit_unpack(input_bytes, m, n, self.ring.bit_unpack_w, packed_len, gamma_2)

    def bit_unpack_z(self, input_bytes, m, n, gamma_1):
        if gamma_1 == (1 << 17):
            packed_len = 576
        elif gamma_1 == (1 << 19):
            packed_len = 640
        else:
            raise ValueError("gamma_1 inválido")
        return self.__bit_unpack(input_bytes, m, n, self.ring.bit_unpack_z, packed_len, gamma_1)

class MatrixDilithium:
    def __init__(self, rows: List[List], ring: PolynomialRingDilithium = None):
        self.m = len(rows)
        self.n = len(rows[0]) if self.m > 0 else 0
        self.data = rows
        self.ring = ring or PolynomialRingDilithium()
        self.parent = self.ring

        self.rows = self.m
        self.cols = self.n

    def __getitem__(self, idx):
        if isinstance(idx, tuple):
            i, j = idx
            return self.data[i][j]
        else:
            return self.data[idx]

    def __setitem__(self, idx, value):
        i, j = idx
        self.data[i][j] = value

    def __len__(self):
        return self.m

    def shape(self):
        return (self.m, self.n)

    def transpose(self):
        transposed = [
            [self.data[i][j] for i in range(self.m)]
            for j in range(self.n)
        ]
        return MatrixDilithium(transposed, self.ring)

    def map(self, f):
        new_data = [[f(self[i, j]) for j in range(self.n)] for i in range(self.m)]
        return MatrixDilithium(new_data, self.ring)

    def to_ntt(self):
        return self.map(lambda x: x.to_ntt())

    def from_ntt(self):
        return self.map(lambda x: x.from_ntt())

    def __add__(self, other):
        if (self.m, self.n) != (other.m, other.n):
            raise ValueError("Dimensões incompatíveis para adição.")
        new_data = [
            [self.data[i][j] + other.data[i][j] for j in range(self.n)]
            for i in range(self.m)
        ]
        return MatrixDilithium(new_data, self.ring)

    def __mul__(self, other):
        if self.n != other.m:
            raise ValueError("Dimensões incompatíveis para multiplicação.")
        result = [[None for _ in range(other.n)] for _ in range(self.m)]
        for i in range(self.m):
            for j in range(other.n):
                acc = self.ring(0)
                for k in range(self.n):
                    acc += self.data[i][k] * other.data[k][j]
                result[i][j] = acc
        return MatrixDilithium(result, self.ring)

    def __repr__(self):
        return "\n".join(["[" + ", ".join(str(p) for p in row) + "]" for row in self.data])
    
    @classmethod
    def from_seed(cls, seed: bytes, rows: int, cols: int, R: PolynomialRingDilithium):

        mat = []
        for i in range(rows):
            row = []
            for j in range(cols):
                poly = R.rejection_sample_ntt_poly(seed, i, j)
                row.append(poly)
            mat.append(row)
        return cls(mat, R)
    
    def sample_vector_ntt(self, seed: bytes, dim: int):

        return MatrixDilithium(
            [[self.ring.rejection_sample_ntt_poly(seed, 0, j)] for j in range(dim)],
            self.ring
        )
    
    def __sub__(self, other):
        if self.m != other.m or self.n != other.n:
            raise ValueError("Dimensões incompatíveis para subtração.")
        result = [
            [self.data[i][j] - other.data[i][j] for j in range(self.n)]
            for i in range(self.m)
        ]
        return MatrixDilithium(result, self.parent)

    
def test_modules():
    print("Teste: Matrizes de polinómios")

    R = PolynomialRingDilithium()
    seed = b"A" * 32

    # A já vem com polinómios em NTT
    A = MatrixDilithium.from_seed(seed, 3, 3, R)

    # Verifica dimensão
    assert A.rows == 3
    assert A.cols == 3

    # Faz from_ntt() diretamente, já que A está em forma NTT
    A_recovered = A.from_ntt()
    for i in range(3):
        for j in range(3):
            assert A[i, j].from_ntt().coeffs == A_recovered[i, j].coeffs

    # Multiplicação A x s
    s = A.sample_vector_ntt(seed, dim=3)  # já vem em NTT
    result = A * s
    assert len(result) == 3
    assert all(len(result[i, 0].coeffs) == 256 for i in range(len(result)))

    print("-> OK")


def run_all_tests():
    test_modules()
    print("Todos os testes passaram com sucesso.")

run_all_tests()