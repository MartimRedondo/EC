import hashlib
from sage.all import *
from typing import Generator
import random

class PolynomialDilithium:
    def __init__(self, parent, coefficients):
        self.parent = parent
        self.q = parent.q
        self.Rq = parent.Rq
        self._coeffs = self._parse_coefficients(coefficients)
        self.poly = self.Rq(self._coeffs)

    @property
    def coeffs(self):
        return self._coeffs

    def _parse_coefficients(self, coeffs):
        n = self.parent.n
        coeffs = list(coeffs)
        if len(coeffs) < n:
            coeffs += [0] * (n - len(coeffs))
        elif len(coeffs) > n:
            coeffs = coeffs[:n]
        return coeffs

    def to_ntt(self):

        coeffs = self.coeffs[:]
        k, l = 0, 128
        zetas = self.parent.ntt_zetas
        while l > 0:
            start = 0
            while start < 256:
                k += 1
                zeta = zetas[k]
                for j in range(start, start + l):
                    t = zeta * coeffs[j + l] % self.q
                    u = coeffs[j]
                    coeffs[j] = (u + t) % self.q
                    coeffs[j + l] = (u - t) % self.q
                start = start + 2 * l
            l >>= 1

        return self.parent(coeffs, is_ntt=True)

    def power_2_round(self, d):
        power_2 = 1 << d
        r1 = []
        r0 = []
        for c in self.coeffs:
            r = int(c) % self.q 
            r0_i = r % power_2  
            r1.append((r - r0_i) >> d)
            r0.append(r0_i)
        return self.parent(r1), self.parent(r0)

    def high_bits(self, alpha):
        def high(c): return ((c + alpha//2) // alpha) % self.q
        coeffs = [high(c) for c in self.coeffs]
        return self.parent(coeffs)

    def low_bits(self, alpha):
        def low(c): return (c - ((c + alpha//2) // alpha) * alpha) % self.q
        coeffs = [low(c) for c in self.coeffs]
        return self.parent(coeffs)

    def decompose(self, alpha):
        high = []
        low = []
        for c in self.coeffs:
            high_c = ((c + alpha//2) // alpha) % self.q
            low_c = (c - high_c * alpha) % self.q
            high.append(high_c)
            low.append(low_c)
        return self.parent(high), self.parent(low)

    def check_norm_bound(self, bound):
        for c in self.coeffs:
            c_int = int(c)
            centered = c_int if c_int < self.q // 2 else c_int - self.q
            if abs(centered) >= bound:
                return True
        return False

    def __add__(self, other):
        if not isinstance(other, PolynomialDilithium):
            raise NotImplementedError("Soma com tipo inválido")
        return self.parent([(a + b) % self.q for a, b in zip(self.coeffs, other.coeffs)])

    def __sub__(self, other):
        if not isinstance(other, PolynomialDilithium):
            raise NotImplementedError("Sub com tipo inválido")
        return self.parent([(a - b) % self.q for a, b in zip(self.coeffs, other.coeffs)])

    def __mul__(self, other):
        if isinstance(other, PolynomialDilithium):
            prod = self.poly * other.poly
            return self.parent(prod.lift().list()[:self.parent.n])
        elif isinstance(other, int):
            return self.parent([(a * other) % self.q for a in self.coeffs])
        else:
            raise TypeError("Multiplicação inválida.")
        
    def __rmul__(self, other):
        return self * other
 
    def __pow__(self, e):
        if not isinstance(e, int):
            raise TypeError("Expoente deve ser inteiro")
        if e < 0:
            raise ValueError("Expoente não pode ser negativo")
        result = self.parent(1)
        for _ in range(e):
            result *= self
        return result

    def __repr__(self):
        terms = []
        for i, c in enumerate(self.coeffs):
            if c == 0:
                continue
            if i == 0:
                terms.append(f"{c}")
            elif i == 1:
                terms.append(f"{'' if c == 1 else c}x")
            else:
                terms.append(f"{'' if c == 1 else c}x^{i}")
        return " + ".join(terms) if terms else "0"
    
    @staticmethod
    def __bit_pack(coeffs, n_bits, n_bytes):

        packed = bytearray()
        buffer = 0
        bits_in_buffer = 0

        for c in coeffs:
            buffer |= (int(c) & ((1 << n_bits) - 1)) << bits_in_buffer
            bits_in_buffer += n_bits
            while bits_in_buffer >= 8:
                packed.append(buffer & 0xFF)
                buffer >>= 8
                bits_in_buffer -= 8

        if bits_in_buffer > 0:
            packed.append(buffer)

        if len(packed) != n_bytes:
            raise ValueError(f"bit_pack falhou: esperado {n_bytes} bytes, obtido {len(packed)}")
        return bytes(packed)

    def bit_pack_t0(self):
        altered = [(1 << 12) - c for c in self.coeffs]
        return self.__bit_pack(altered, 13, 416)

    def bit_pack_t1(self):
        for c in self.coeffs:
            if not (0 <= c < (1 << 10)):
                raise ValueError(f"Coeficiente inválido para t1: {c}")

        return self.__bit_pack(self.coeffs, 10, 320)


    def bit_pack_s(self, eta):
        altered = [(eta - c) % self.q for c in self.coeffs]
        if eta == 2:
            return self.__bit_pack(altered, 3, 96)
        elif eta == 4:
            return self.__bit_pack(altered, 4, 128)
        else:
            raise ValueError("eta inválido")

    def bit_pack_w(self, gamma_2):
        if gamma_2 == 95232:
            return self.__bit_pack(self.coeffs, 6, 192)
        elif gamma_2 == 261888:
            return self.__bit_pack(self.coeffs, 4, 128)
        else:
            raise ValueError("gamma_2 inválido")

    def bit_pack_z(self, gamma_1):
        altered = [(gamma_1 - c) % self.q for c in self.coeffs]
        if gamma_1 == (1 << 17):
            return self.__bit_pack(altered, 18, 576)
        elif gamma_1 == (1 << 19):
            return self.__bit_pack(altered, 20, 640)
        else:
            raise ValueError("gamma_1 inválido")
        
    def make_hint(self, other, alpha):
        """
        h_i = 1 se high_bits(r_i + alpha // 2) != high_bits(z_i + alpha // 2), senão 0
        """
        def high_bits(c): return ((c + alpha//2) // alpha) % self.q

        hints = [
            1 if high_bits(r) != high_bits(z) else 0
            for r, z in zip(self.coeffs, other.coeffs)
        ]
        return self.parent(hints)

    def make_hint_optimised(self, other, alpha):
        """
        Optimized version of make_hint:
        h_i = 1 se low_bits(z_i) > alpha / 2 ou se high_bits(r_i) != high_bits(z_i)
        """
        def high_bits(c): return ((c + alpha//2) // alpha) % self.q
        def low_bits(c): return (c - ((c + alpha//2) // alpha) * alpha) % self.q

        hints = [
            1 if low_bits(z) > (alpha // 2) or high_bits(r) != high_bits(z) else 0
            for r, z in zip(self.coeffs, other.coeffs)
        ]
        return self.parent(hints)

    def use_hint(self, hint_poly, alpha):
        """
        h_i == 0 → return high_bits(self_i)
        h_i == 1 → recomputar usando uso do hint
        """
        def high_bits(c): return ((c + alpha//2) // alpha) % self.q
        def low_bits(c): return (c - ((c + alpha//2) // alpha) * alpha) % self.q

        result = []
        for h, r in zip(hint_poly.coeffs, self.coeffs):
            if h == 0:
                result.append(high_bits(r))
            else:
                r_high = high_bits(r)
                r_low = low_bits(r)
                if r_low > alpha // 2:
                    r_high = (r_high + 1) % self.q
                result.append(r_high)
        return self.parent(result)
    
    def from_ntt(self):
        # Se não está em NTT, retorna a si próprio
        return self
    
    def is_zero(self):
        return all(c == 0 for c in self.coeffs)

    def is_constant(self):
        return all(c == 0 for c in self.coeffs[1:])

    def reduce_coefficients(self):
        reduced = [c % self.q for c in self.coeffs]
        return self.parent(reduced)

    def __eq__(self, other):
        if not isinstance(other, PolynomialDilithium):
            return False
        return self.coeffs == other.coeffs

    def __neg__(self):
        negated = [(-c) % self.q for c in self.coeffs]
        return self.parent(negated)

    def __getitem__(self, i):
        return self.coeffs[i]

    @property
    def is_ntt(self):
        return isinstance(self, PolynomialDilithiumNTT)

    
class PolynomialDilithiumNTT(PolynomialDilithium):
    def __init__(self, parent, coefficients):
        self.parent = parent
        self.q = parent.q
        self.Rq = parent.Rq
        self._coeffs = list(coefficients)
        self.poly = self.Rq(coefficients)

    def to_ntt(self):
        raise TypeError("Este polinómio já está em forma NTT")

    def from_ntt(self):

        coeffs = self.coeffs[:]
        l, k = 1, 256
        zetas = self.parent.ntt_zetas
        while l < 256:
            start = 0
            while start < 256:
                k -= 1
                zeta = -zetas[k]
                for j in range(start, start + l):
                    t = coeffs[j]
                    u = coeffs[j + l]
                    coeffs[j] = (t + u) % self.q
                    coeffs[j + l] = ((t - u) * zeta) % self.q
                start = start + 2 * l
            l <<= 1

        f = self.parent.ntt_f
        coeffs = [(c * f) % self.q for c in coeffs]
        return self.parent(coeffs, is_ntt=False)

    def __add__(self, other):
        if not isinstance(other, PolynomialDilithium):
            raise NotImplementedError("Soma com tipo inválido")
        return self.parent([(a + b) % self.q for a, b in zip(self.coeffs, other.coeffs)], is_ntt=True)

    def __sub__(self, other):
        if not isinstance(other, PolynomialDilithium):
            raise NotImplementedError("Sub com tipo inválido")
        return self.parent([(a - b) % self.q for a, b in zip(self.coeffs, other.coeffs)], is_ntt=True)

    def __mul__(self, other):
        if isinstance(other, PolynomialDilithium):
            if isinstance(other, PolynomialDilithiumNTT):
                return self.parent([(a * b) % self.q for a, b in zip(self.coeffs, other.coeffs)], is_ntt=True)
            else:
                raise TypeError("Multiplicação entre NTT e não-NTT não suportada")
        elif isinstance(other, int):
            return self.parent([(c * other) % self.q for c in self.coeffs], is_ntt=True)
        else:
            raise TypeError("Multiplicação inválida em NTT")

class PolynomialRingDilithium:
    def __init__(self):
        self.q = 8380417
        self.n = 256
        self.Zq = Integers(self.q)
        self.R = PolynomialRing(self.Zq, 'x')  # R = Z_q[x]
        self.x = self.R.gen()
        self.Rq = self.R.quotient(self.x**self.n + 1, 'x')
        self.element = PolynomialDilithium
        self.element_ntt = PolynomialDilithiumNTT

        root_of_unity = 1753
        self.ntt_zetas = [
            power_mod(root_of_unity, self._bit_reverse(i, 8), self.q)
            for i in range(self.n)
        ]
        self.ntt_f = inverse_mod(self.n, self.q)

    def _bit_reverse(self, i, k):
        return int(bin(i)[2:].zfill(k)[::-1], 2)

    def __call__(self, coeffs, is_ntt=False):
        if isinstance(coeffs, int):
            coeffs = [coeffs]
        if not isinstance(coeffs, list):
            raise TypeError("Esperada lista de coeficientes")
        if len(coeffs) > self.n:
            raise ValueError("Lista demasiado longa")
        element = self.element_ntt if is_ntt else self.element
        return element(self, coeffs)

    def __bit_unpack(self, input_bytes, n_bits):
        total_bits = self.n * n_bits
        if len(input_bytes) * 8 != total_bits:
            raise ValueError("Número de bytes não corresponde ao número de bits esperado")
        r = int.from_bytes(input_bytes, "little")
        mask = (1 << n_bits) - 1
        return [(r >> (n_bits * i)) & mask for i in range(self.n)]

    def bit_unpack_t0(self, input_bytes):
        altered = self.__bit_unpack(input_bytes, 13)
        coeffs = [(1 << 12) - c for c in altered]
        return self(coeffs)

    def bit_unpack_t1(self, packed_bytes):
        assert len(packed_bytes) == 320, f"Esperado 320 bytes, recebido {len(packed_bytes)}"

        coeffs = []
        acc = 0
        bits = 0
        for b in packed_bytes:
            acc |= b << bits
            bits += 8
            while bits >= 10:
                coeffs.append(acc & 0x3FF)
                acc >>= 10
                bits -= 10

        if len(coeffs) != self.n:
            raise ValueError(f"Número incorreto de coeficientes: {len(coeffs)} (esperado {self.n})")

        return self(coeffs)

    def bit_unpack_s(self, input_bytes, eta):
        if eta == 2:
            altered = self.__bit_unpack(input_bytes, 3)
        elif eta == 4:
            altered = self.__bit_unpack(input_bytes, 4)
        else:
            raise ValueError("eta inválido")
        coeffs = [(eta - c) % self.q for c in altered]
        return self(coeffs)

    def bit_unpack_w(self, input_bytes, gamma_2):
        if gamma_2 == 95232:
            bits = 6
        elif gamma_2 == 261888:
            bits = 4
        else:
            raise ValueError("gamma_2 inválido")
        coeffs = self.__bit_unpack(input_bytes, bits)
        return self(coeffs)

    def bit_unpack_z(self, input_bytes, gamma_1):
        if gamma_1 == (1 << 17):
            bits = 18
        elif gamma_1 == (1 << 19):
            bits = 20
        else:
            raise ValueError("gamma_1 inválido")
        altered = self.__bit_unpack(input_bytes, bits)
        coeffs = [(gamma_1 - c) % self.q for c in altered]
        return self(coeffs)

    def _shake256_stream(self, seed: bytes) -> Generator[int, None, None]:
        xof = hashlib.shake_256()
        xof.update(seed)
        while True:
            for b in xof.digest(4096):  # Yields one byte at a time
                yield b

    def _shake128_stream(self, seed: bytes) -> Generator[int, None, None]:
        xof = hashlib.shake_128()
        xof.update(seed)
        while True:
            for b in xof.digest(4096):
                yield b

    def sample_in_ball(self, seed: bytes, tau: int):
        xof = self._shake256_stream(seed)
        sign_bytes = bytes([next(xof) for _ in range(8)])
        sign_int = int.from_bytes(sign_bytes, "little")

        coeffs = [0] * self.n
        selected = set()

        while len(selected) < tau:
            j = next(xof)
            if j < self.n and j not in selected:
                selected.add(j)

        selected = list(selected)
        for idx in selected:
            coeffs[idx] = 1 - 2 * (sign_int & 1)
            sign_int >>= 1

        return self(coeffs)

    def rejection_sample_ntt_poly(self, rho: bytes, i: int, j: int):
        seed = rho + bytes([i, j])
        xof = self._shake128_stream(seed)
        coeffs = []
        while len(coeffs) < self.n:
            b = bytes([next(xof) for _ in range(3)])
            val = int.from_bytes(b, "little") & 0x7FFFFF
            if val < self.q:
                coeffs.append(val)
        return self(coeffs, is_ntt=True)

    def rejection_bounded_poly(self, rho_prime: bytes, i: int, eta: int):
        def decode_half_byte(j, eta):
            if eta == 2 and j < 15:
                return 2 - (j % 5)
            elif eta == 4 and j < 9:
                return 4 - j
            return None

        seed = rho_prime + int(i).to_bytes(2, "little")
        xof = self._shake256_stream(seed)
        coeffs = []
        while len(coeffs) < self.n:
            j = next(xof)
            c0 = decode_half_byte(j & 0x0F, eta)
            if c0 is not None:
                coeffs.append(c0)
            if len(coeffs) < self.n:
                c1 = decode_half_byte(j >> 4, eta)
                if c1 is not None:
                    coeffs.append(c1)
        return self(coeffs)

    def sample_mask_polynomial(self, rho_prime: bytes, i: int, kappa: int, gamma_1: int):

        bit_count = 18 if gamma_1 == (1 << 17) else 20
        total_bytes = (self.n * bit_count) // 8

        seed = rho_prime + int(kappa + i).to_bytes(2, "little")
        xof = hashlib.shake_256(seed).digest(total_bytes)
        r = int.from_bytes(xof, "little")
        mask = (1 << bit_count) - 1
        coeffs = [gamma_1 - ((r >> (bit_count * i)) & mask) for i in range(self.n)]

        return self(coeffs)
    
    def sample_eta(self, seed: bytes, nonce: int, eta: int = 2):
        shake = hashlib.shake_256()
        shake.update(seed + nonce.to_bytes(2, "little"))
        buf = bytearray(shake.digest(512))  # podemos ajustar o tamanho conforme necessário

        coeffs = []
        i = 0
        while len(coeffs) < self.n and i < len(buf):
            byte = buf[i]
            if eta == 2 and byte < 243:
                for j in range(0, 8, 4):
                    val = (byte >> j) & 0x03
                    coeffs.append(val - 2)
                    if len(coeffs) == self.n:
                        break
            elif eta == 4 and byte < 248:
                for j in range(0, 8, 2):
                    val = (byte >> j) & 0x03
                    coeffs.append(val - 2)
                    if len(coeffs) == self.n:
                        break
            i += 1

        coeffs = coeffs[:self.n]  # cortar se passou do limite
        return PolynomialDilithium(self, coeffs)
    
    def gen(self):
        return self([0, 1])

    def random_element(self):
        import random
        coeffs = [ZZ.random_element(0, self.q) for _ in range(self.n)]
        return self(coeffs)

    def __repr__(self):
        return f"Univariate Polynomial Ring in x over Z/{self.q}Z (modulus x^{self.n} + 1)"

    def random(self):
        import random
        return self([random.randint(0, self.q - 1) for _ in range(self.n)])




def test_polynomial_ring():
    print("Teste: Criação do anel")
    R = PolynomialRingDilithium()
    assert R.q == 8380417
    assert R.n == 256
    print("-> OK")

def test_arithmetic():
    print("Teste: Aritmética básica")
    R = PolynomialRingDilithium()
    a = R([1]*256)
    b = R([2]*256)
    c = a + b
    d = b - a
    e = a * 3
    assert all(x == 3 for x in c.coeffs)
    assert all(x == 1 for x in d.coeffs)
    assert all(x == 3 for x in e.coeffs)
    print("-> OK")

def test_ntt_roundtrip():
    print("Teste: NTT -> inversa")
    R = PolynomialRingDilithium()
    a = R([i for i in range(256)])
    a_ntt = a.to_ntt()
    a_recovered = a_ntt.from_ntt()
    diff = [(x - y) % R.q for x, y in zip(a.coeffs, a_recovered.coeffs)]
    assert sum(diff) == 0
    print("-> OK")

def test_bit_packing():
    print("Teste: bit_pack / unpack")
    R = PolynomialRingDilithium()
    poly = R([i % 4096 for i in range(256)])
    packed = poly.bit_pack_t0()
    unpacked = R.bit_unpack_t0(packed)
    diff = [(x - y) % R.q for x, y in zip(poly.coeffs, unpacked.coeffs)]
    assert sum(diff) == 0
    print("-> OK")

def test_sampling():
    print("Teste: sample_in_ball / rejection")
    R = PolynomialRingDilithium()
    seed = b"A" * 32
    p1 = R.sample_in_ball(seed, tau=39)
    print(f"sample_in_ball -> len: {len(p1.coeffs)}, nonzero: {sum(1 for c in p1.coeffs if c != 0)}")
    assert len(p1.coeffs) == 256
    assert sum(abs(c) for c in p1.coeffs) == 39

    p2 = R.rejection_bounded_poly(seed, 0, eta=2)
    assert len(p2.coeffs) == 256
    assert all(-2 <= c <= 2 for c in p2.coeffs)

    p3 = R.rejection_sample_ntt_poly(seed, 0, 0)
    assert len(p3.coeffs) == 256
    assert all(0 <= c < R.q for c in p3.coeffs)

    print("-> OK")

def run_all_tests():
    test_polynomial_ring()
    test_arithmetic()
    test_ntt_roundtrip()
    test_bit_packing()
    test_sampling()
    print("✅ Todos os testes passaram com sucesso.")


if __name__ == "__main__":
    from polynomials_sage import PolynomialRingDilithium
    from random import randint

    R = PolynomialRingDilithium()

    # Gera coeficientes aleatórios de 10 bits (entre 0 e 1023)
    coeffs = [randint(0, (1 << 10) - 1) for _ in range(R.n)]
    p = R(coeffs)

    # Faz o empacotamento e descompactamento
    packed = p.bit_pack_t1()
    unpacked = R.bit_unpack_t1(packed)

    print("bit_pack_t1 e bit_unpack_t1 executados.")
    diff = [(a - b) % R.q for a, b in zip(p.coeffs, unpacked.coeffs)]
    diffs = sum(1 for d in diff if d != 0)
    print(f"Número de coeficientes diferentes: {diffs}")
    assert all(d == 0 for d in diff), "Erro: bit_pack_t1() e bit_unpack_t1() não são inversos perfeitos"

    print("Teste de t1 packing/unpacking passou com sucesso!")


