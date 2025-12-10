import unittest
from dilithium_sage import Dilithium

class TestDilithiumSage(unittest.TestCase):
    def setUp(self):
        self.d = Dilithium()
        self.pk, self.sk = self.d.keygen()
        self.msg = b"mensagem de teste"

    def test_keygen_signature_verification(self):
        sig = self.d.sign(self.msg, self.sk, self.pk)
        self.assertTrue(self.d.verify(self.msg, sig, self.pk), "A assinatura devia ser válida")

    def test_altered_message_invalid_signature(self):
        sig = self.d.sign(self.msg, self.sk, self.pk)
        altered_msg = self.msg + b"!"
        self.assertFalse(self.d.verify(altered_msg, sig, self.pk), "A assinatura não devia ser válida para mensagem alterada")

    def test_corrupted_signature(self):
        sig = bytearray(self.d.sign(self.msg, self.sk, self.pk))
        sig[10] ^= 0xFF  # Inverter alguns bits
        self.assertFalse(self.d.verify(self.msg, bytes(sig), self.pk), "A assinatura não devia ser válida se for corrompida")

    def test_multiple_signatures(self):
        sigs = [self.d.sign(self.msg, self.sk, self.pk) for _ in range(5)]
        for sig in sigs:
            self.assertTrue(self.d.verify(self.msg, sig, self.pk), "Assinatura devia ser válida")

    def test_signature_length(self):
        sig = self.d.sign(self.msg, self.sk, self.pk)
        self.assertEqual(len(sig), self.d.siglen, "Comprimento da assinatura deve ser igual a siglen")

    def test_signature_uniqueness(self):
        sig1 = self.d.sign(self.msg, self.sk, self.pk)
        sig2 = self.d.sign(self.msg, self.sk, self.pk)
        self.assertNotEqual(sig1, sig2, "Assinaturas devem ser distintas (aleatoriedade)")

if __name__ == "__main__":
    unittest.main()
