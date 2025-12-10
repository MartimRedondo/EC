from sage.all import *
import hashlib
import os

from polynomials_sage import PolynomialRingDilithium, PolynomialDilithium
from modules_sage import MatrixDilithium

CRYPTO_BYTES = 2416
PUBLIC_KEY_BYTES = 1312
SECRET_KEY_BYTES = 2528

class Dilithium:
    def __init__(self, mode=2):
        if mode == 2:
            self.k = 4
            self.l = 4
            self.eta = 2
            self.gamma_1 = 2**17
            self.gamma_2 = self.gamma_1 // 2
            self.beta = 78
            self.omega = 80
            self.tau = 39
        elif mode == 3:
            self.k = 6
            self.l = 5
            self.eta = 4
            self.gamma_1 = 2**19
            self.gamma_2 = self.gamma_1 // 4
            self.beta = 196
            self.omega = 55
            self.tau = 49
        else:
            raise ValueError("Modo inválido para Dilithium")

        self.seed_bytes = 32
        self.crh_bytes = 48

        self.siglen = CRYPTO_BYTES
        self.R = PolynomialRingDilithium()

    def _shake256(self, data: bytes, outlen: int) -> bytes:
        return hashlib.shake_256(data).digest(outlen)

    def _shake128(self, data: bytes, outlen: int) -> bytes:
        return hashlib.shake_128(data).digest(outlen)

    def _expand_matrix_from_seed(self, seed_A: bytes) -> MatrixDilithium:
        return MatrixDilithium.from_seed(seed_A, self.k, self.l, self.R)

    def _expand_vector_from_seed(self, seed: bytes) -> MatrixDilithium:
        return MatrixDilithium(
            [[self.R.sample_eta(seed, i, self.eta)] for i in range(self.l)], self.R
        )

    def _pack_pk(self, rho: bytes, t1_matrix: MatrixDilithium) -> bytes:
        packed = bytearray(rho)
        for i in range(self.k):
            packed.extend(t1_matrix[i, 0].bit_pack_t1())
        return bytes(packed)

    def _unpack_pk(self, pk: bytes):
        rho = pk[:self.seed_bytes]
        t1_bytes = pk[self.seed_bytes:]
        t1 = [
            [self.R.bit_unpack_t1(t1_bytes[i * 320 : (i + 1) * 320])]
            for i in range(self.k)
        ]
        return rho, MatrixDilithium(t1, self.R)

    def _pack_sk(self, rho, key, tr, s1, s2, t0):
        sk = bytearray()
        sk.extend(rho)
        sk.extend(key)
        sk.extend(tr)
        for i in range(self.l):
            sk.extend(s1[i, 0].bit_pack_s(self.eta))
        for i in range(self.k):
            sk.extend(s2[i, 0].bit_pack_s(self.eta))
        for i in range(self.k):
            sk.extend(t0[i, 0].bit_pack_t0())
        return bytes(sk)

    def _generate_challenge(self, c_bytes: bytes):
        n = self.R.n
        c = [0] * n
        signs = int.from_bytes(c_bytes[:8], "little")
        pos = 8
        cnt = 0
        print("Challenge generation started")

        while cnt < self.tau:
            if pos >= len(c_bytes):
                c_bytes += self._shake256(c_bytes, 64)
            b = c_bytes[pos] % n
            pos += 1
            if c[b] != 0:
                continue
            c[b] = 1 if (signs & 1) == 0 else -1
            signs >>= 1
            cnt += 1
            if cnt % 10 == 0:
                print(f"{cnt}/{self.tau} non-zero coefficients set")

        print("Challenge generated")
        return self.R(c)

    def keygen(self):
        seed = os.urandom(self.seed_bytes)
        buf = self._shake256(seed, 2 * self.seed_bytes + self.crh_bytes)
        rho = buf[0:self.seed_bytes]
        rho_prime = buf[self.seed_bytes:self.seed_bytes + self.crh_bytes]
        key = buf[self.seed_bytes + self.crh_bytes:]

        A = self._expand_matrix_from_seed(rho)
        s1 = self._expand_vector_from_seed(rho_prime)
        s2 = MatrixDilithium(
            [[self.R.sample_eta(rho_prime, i + self.l, self.eta)] for i in range(self.k)],
            self.R
        )

        s1_ntt = s1.to_ntt()
        t = A * s1_ntt
        t = t.from_ntt() + s2

        t0 = t.map(lambda p: p.power_2_round(13)[1])  # parte baixa
        t1 = t.map(lambda p: p.power_2_round(13)[0])  # parte alta (deve ir para a chave pública)

        pk = self._pack_pk(rho, t1)
        tr = self._shake256(pk, self.crh_bytes)
        sk = self._pack_sk(rho, key, tr, s1, s2, t0)

        return pk, sk

    def sign(self, message: bytes, sk: bytes, pk: bytes) -> bytes:
        rho = sk[:self.seed_bytes]
        print("RHO - SIGN:")
        print(rho.hex())
        # Extração de t1 da chave pública
        t1 = MatrixDilithium(
            [[self.R.bit_unpack_t1(pk[self.seed_bytes + i * 320 : self.seed_bytes + (i + 1) * 320])] for i in range(self.k)],
            self.R
        )
        offset = self.seed_bytes * 2 + self.crh_bytes

        s1 = MatrixDilithium(
            [[self.R.bit_unpack_s(sk[offset + i * 96 : offset + (i + 1) * 96], self.eta)] for i in range(self.l)], self.R
        )
        offset += self.l * 96

        s2 = MatrixDilithium(
            [[self.R.bit_unpack_s(sk[offset + i * 96 : offset + (i + 1) * 96], self.eta)] for i in range(self.k)], self.R
        )
        offset += self.k * 96

        t0 = MatrixDilithium(
            [[self.R.bit_unpack_t0(sk[offset + i * 416 : offset + (i + 1) * 416])] for i in range(self.k)], self.R
        )

        mu = self._shake256(self._shake256(pk, self.crh_bytes) + message, self.crh_bytes)
        A = self._expand_matrix_from_seed(rho)

        attempt = 0
        while True:
            attempt += 1
            print(f"Tentativa de assinatura {attempt}")

            y = MatrixDilithium(
                [[self.R.sample_eta(os.urandom(32), 100 + i, self.eta)] for i in range(self.l)], self.R
            ).to_ntt()

            w = A * y
            w1 = w.map(lambda p: p.power_2_round(9)[1])

            print("w1 assinatura exemplo:", [p.coeffs[:4] for row in w1.data for p in row][:2])

            w1_bytes = b"".join([p.bit_pack_t1() for row in w1.data for p in row])
            c_hash = self._shake256(mu + w1_bytes, self.seed_bytes)

            c_poly = self._generate_challenge(c_hash)
            z = y.from_ntt() + s1.map(lambda p: c_poly * p)

            if all(not p.check_norm_bound(self.gamma_1 - self.beta) for row in z.data for p in row):
                print("Assinatura aceite")
                break
            else:
                print("Assinatura rejeitada")

        z_bytes = b"".join([p.bit_pack_z(self.gamma_1) for row in z.data for p in row])
        c_bytes = c_hash
        h_bytes = b"\x00" * self.omega

        self._debug_last_sign = {
            "mu": mu,
            "w1_bytes": w1_bytes,
            "c_hash": c_hash,
            "w1": w1,  # guarda a matriz completa
            "A": [[p.coeffs[:4] for p in row] for row in A.data]    # guarda a matriz completa
        }

        self.pk_cached = pk
        self.t1_cached = t1

        return c_bytes + z_bytes + h_bytes

    def verify(self, message: bytes, signature: bytes, pk: bytes) -> bool:
        print("Início da verificação")
        c_bytes = signature[:self.seed_bytes]
        print("c_bytes extraído:", c_bytes.hex())
        z_bytes = signature[self.seed_bytes:self.seed_bytes + self.l * 576]
        print("z_bytes extraído com comprimento:", len(z_bytes))

        z = MatrixDilithium(
            [[self.R.bit_unpack_z(z_bytes[i * 576:(i + 1) * 576], self.gamma_1)] for i in range(self.l)], self.R
        )
        print("z descompactado")

        rho, t1 = self._unpack_pk(pk)
        print("RHO - VERIFY")
        print(rho.hex())
        print("Chave pública descompactada")
        mu = self._shake256(self._shake256(pk, self.crh_bytes) + message, self.crh_bytes)
        print("Digest mu computado")

        debug = getattr(self, "_debug_last_sign", None)
        A = self._expand_matrix_from_seed(rho)
        print("Matriz A expandida")

        if debug and "A" in debug:
            A_coeffs = [[p.coeffs[:4] for p in row] for row in A.data]
            same = A_coeffs == debug["A"]
            print("Matriz A igual à usada na assinatura?", same)
            if not same:
                for i in range(min(len(A_coeffs), len(debug["A"]))):
                    for j in range(min(len(A_coeffs[i]), len(debug["A"][i]))):
                        if A_coeffs[i][j] != debug["A"][i][j]:
                            print(f"Diferença em A[{i}][{j}]: {debug['A'][i][j]} vs {A_coeffs[i][j]}")
                            break

        if debug:
            print("Comparação com valores da assinatura:")
            print("mu assinada == mu verificada?", debug["mu"] == mu)

        c_poly = self._generate_challenge(c_bytes)
        print("Challenge c_poly gerado")

        ct1 = t1.map(lambda p: c_poly * p)
        z_ntt = z.to_ntt()
        w_prime = A * z_ntt - ct1
        print("w_prime computado")

        w1_prime = self._reconstruct_w1(z, c_poly, A, t1)
        print("w1_prime verificação exemplo:", [p.coeffs[:4] for row in w1_prime.data for p in row][:2])

        if debug and "w1" in debug:
            flat_w1 = [p.coeffs[:4] for row in w1_prime.data for p in row]
            for i in range(min(len(debug["w1"]), len(flat_w1))):
                if flat_w1[i] != debug["w1"][i]:
                    print(f"Diferença em w1[{i}]: esperado {debug['w1'][i]}, obtido {flat_w1[i]}")
                    break

        w1_bytes = b"".join([p.bit_pack_t1() for row in w1_prime.data for p in row])
        print("w1_prime empacotado")

        if debug:
            matches = debug["w1_bytes"] == w1_bytes
            print("w1_bytes assinada == recomputada?", matches)
            if not matches:
                mismatches = sum(a != b for a, b in zip(debug["w1_bytes"], w1_bytes))
                print(f"w1_bytes tem {mismatches} bytes diferentes")

        c_prime = self._shake256(mu + w1_bytes, self.seed_bytes)
        print("c_prime:", c_prime.hex())
        print("c_bytes:", c_bytes.hex())

        result = c_prime == c_bytes
        print("Resultado da verificação:", "VÁLIDA" if result else "INVÁLIDA")

        if debug:
            from difflib import SequenceMatcher
            matcher = SequenceMatcher(None, debug["w1_bytes"], w1_bytes)
            print(f"Percentual de match entre w1_bytes: {matcher.ratio()*100:.2f}%")

        return result


    def _reconstruct_w1(self, z, c, A, t1):
        """
        Reconstrói w1 = HighBits(A * z - c * t1)
        """
        Az = A * z.to_ntt()             # Az in NTT
        ct1 = t1.map(lambda p: c * p)   # c * t1
        w = Az - ct1                    # A*z - c*t1
        w = w.from_ntt()                # Volta para domínio normal
        w1 = w.map(lambda p: p.power_2_round(9)[1])  # Extrai high bits (t1)
        return w1


if __name__ == "__main__":
    dilithium = Dilithium(mode=2)
    pk, sk = dilithium.keygen()
    print(f"Chave pública gerada com {len(pk)} bytes.")
    print(f"Chave privada gerada com {len(sk)} bytes.")

    mensagem = b"mensagem de teste"
    assinatura = dilithium.sign(mensagem, sk, pk)
    print(f"Assinatura gerada com {len(assinatura)} bytes.")

    resultado = dilithium.verify(mensagem, assinatura, pk)
    print(f"Assinatura válida? {resultado}")

