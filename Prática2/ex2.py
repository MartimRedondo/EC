import hashlib, os 
from pickle import dumps
from sage.all import *

import time


class Ed(object):
    def __init__(self,p, a, d , ed = None):
        assert a != d and is_prime(p) and p > 3
        K        = GF(p) 
  
        A =  2*(a + d)/(a - d)
        B =  4/(a - d)
    
        alfa = A/(3*B) ; s = B

        a4 =  s**(-2) - 3*alfa**2
        a6 =  -alfa**3 - a4*alfa
        
        self.K = K
        self.constants = {'a': a , 'd': d , 'A':A , 'B':B , 'alfa':alfa , 's':s , 'a4':a4 , 'a6':a6 }
        self.EC = EllipticCurve(K,[a4,a6]) 
        
        if ed != None:
            self.L = ed['L']
            self.P = self.ed2ec(ed['Px'],ed['Py'])  # gerador do gru
        else:
            self.gen()
    
    def order(self):
        # A ordem prima "n" do maior subgrupo da curva, e o respetivo cofator "h" 
        oo = self.EC.order()
        n,_ = list(factor(oo))[-1]
        return (n,oo//n)
    
    
    def gen(self):
        L, h = self.order()       
        P = O = self.EC(0)
        while L*P == O:
            P = self.EC.random_element()
        self.P = h*P ; self.L = L

  
    def is_edwards(self, x, y):
        a = self.constants['a'] ; d = self.constants['d']
        x2 = x**2 ; y2 = y**2
        return a*x2 + y2 == 1 + d*x2*y2

    def ed2ec(self,x,y):      ## mapeia Ed --> EC
        if (x,y) == (0,1):
            return self.EC(0)
        z = (1+y)/(1-y) ; w = z/x
        alfa = self.constants['alfa']; s = self.constants['s']
        return self.EC(z/s + alfa , w/s)
    
    def ec2ed(self,P):        ## mapeia EC --> Ed
        if P == self.EC(0):
            return (0,1)
        x,y = P.xy()
        alfa = self.constants['alfa']; s = self.constants['s']
        u = s*(x - alfa) ; v = s*y
        return (u/v , (u-1)/(u+1))


class ed(object):
    def __init__(self,pt=None,curve=None,x=None,y=None):
        if pt != None:
            self.curve = pt.curve
            self.x = pt.x ; self.y = pt.y ; self.w = pt.w
        else:
            assert isinstance(curve,Ed) and curve.is_edwards(x,y)
            self.curve = curve
            self.x = x ; self.y = y ; self.w = x*y
    
    def eq(self,other):
        return self.x == other.x and self.y == other.y
    
    def copy(self):
        return ed(curve=self.curve, x=self.x, y=self.y)
    
    def zero(self):
        return ed(curve=self.curve,x=0,y=1)
    
    def sim(self):
        return ed(curve=self.curve, x= -self.x, y= self.y)
    
    def soma(self, other):
        a = self.curve.constants['a']; d = self.curve.constants['d']
        delta = d*self.w*other.w
        self.x, self.y  = (self.x*other.y + self.y*other.x)/(1+delta), (self.y*other.y - a*self.x*other.x)/(1-delta)
        self.w = self.x*self.y
        
    def duplica(self):
        a = self.curve.constants['a']; d = self.curve.constants['d']
        delta = d*(self.w)**2
        self.x, self.y = (2*self.w)/(1+delta) , (self.y**2 - a*self.x**2)/(1 - delta)
        self.w = self.x*self.y
        
    def mult(self, n):
        m = Mod(n,self.curve.L).lift().digits(2)   ## obter a representação binária do argumento "n"
        Q = self.copy() ; A = self.zero()
        for b in m:
            if b == 1:
                A.soma(Q)
            Q.duplica()
        return A
            




# ECDSA Implementation following FIPS 186-5 standard
class ECDSA:
    def __init__(self, curve_name):
        # Define curves parameters
        if curve_name == "ed25519":

            p = 2**255-19
            K = GF(p)
            a = K(-1)
            d = -K(121665)/K(121666)

            ed25519 = {
                'b'  : 256,
                'Px' : K(15112221349535400772501151409588531511454012693041857206046113283949847762202),
                'Py' : K(46316835694926478169428394003475163141307993866256225615783033603165251855960),
                'L'  : ZZ(2**252 + 27742317777372353535851937790883648493), ## ordem do subgrupo primo
                'n'  : 254,
                'h'  : 8
            }     

            self.E = Ed(p,a,d,ed=ed25519)
            self.G = ed(curve=self.E, x=ed25519['Px'], y=ed25519['Py'])
            self.b = ed25519['b']
            self.n = ed25519['L']
            self.p = p
            self.c = 3
            self.securitybits = 128
            
        elif curve_name == "ed448":

            p = 2**448 - 2**224 - 1
            K = GF(p)
            a = K(1)
            d = K(-39081)

            ed448= {
            'b'  : 456,     ## tamanho das assinaturas e das chaves públicas
            'Px' : K(224580040295924300187604334099896036246789641632564134246125461686950415467406032909029192869357953282578032075146446173674602635247710) ,
            'Py' : K(298819210078481492676017930443930673437544040154080242095928241372331506189835876003536878655418784733982303233503462500531545062832660) ,                                          
            'L'  : ZZ(2**446 - 13818066809895115352007386748515426880336692474882178609894547503885) ,
            'n'  : 447,     ## tamanho dos segredos: os dois primeiros bits são 0 e o último é 1.
            'h'  : 4        ## cofactor
            }

            self.E = Ed(p,a,d,ed=ed448)
            self.G = ed(curve=self.E, x=ed448['Px'], y=ed448['Py'])
            self.b = ed448['b']
            self.n = ed448['L']
            self.p = p
            self.c = 2
            self.securitybits = 224

        else:
            raise ValueError("Curve name must be 'ed25519' or 'ed448'")

        self.curve = curve_name


    # hash function for each curve ED2556 and ED448
    def hash(self,data):
        if self.curve == 'ed25519':
            return hashlib.sha512(data).digest()
        else:
            return hashlib.shake_256(data).digest(912//8)


    def encode_point(self, P):
        x_int = int(P.x)
        y_int = int(P.y)

        if self.curve == "ed25519":
            y_bytes = y_int.to_bytes(32, byteorder='little')

            # Ensure MSB of last byte is 0
            y_bytes = bytearray(y_bytes)
            y_bytes[31] &= 0x7F  # Clear the MSB

            # Set the MSB to the LSB of x
            if x_int & 1:
                y_bytes[31] |= 0x80

            return bytes(y_bytes)

        else:  # ed448
            y_bytes = y_int.to_bytes(57, byteorder='little')

            # Ensure last byte is 0
            y_bytes[56] = 0

            # Set the MSB of second-to-last byte to LSB of x
            if x_int & 1:
                y_bytes[55] |= 0x80

            return bytes(y_bytes)


    def tonelli_shanks(self, n, p):
        """Tonelli-Shanks algorithm to find the square root of n mod p."""
        if n == 0:
            return 0
        if p == 2:
            return n % 2

        # Check if n is a quadratic residue mod p
        if legendre_symbol(n, p) != 1:
            raise ValueError(f"{n} is not a quadratic residue modulo {p}")

        # If p ≡ 3 (mod 4), we can use a faster method
        if p % 4 == 3:
            return power_mod(n, (p + 1) // 4, p)

        # Tonelli-Shanks algorithm for general primes
        s = 0
        q = p - 1
        while q % 2 == 0:
            s += 1
            q //= 2

        z = 2
        while legendre_symbol(z, p) != -1:
            z += 1

        m = s
        c = power_mod(z, q, p)
        t = power_mod(n, q, p)
        r = power_mod(n, (q + 1) // 2, p)

        while t != 0 and t != 1:
            t2i = t
            i = 0
            for i in range(1, m):
                t2i = power_mod(t2i, 2, p)
                if t2i == 1:
                    break
            if i == m:
                raise ValueError("Tonelli-Shanks failed")

            b = power_mod(c, 2**(m - i - 1), p)
            m = i
            c = power_mod(b, 2, p)
            t = (t * b * b) % p
            r = (r * b) % p

        return r


    def decode_point(self, encoded):
        a = self.E.constants['a']
        d = self.E.constants['d']
        p = self.p

        # Interpret octet string as integer
        integer = int.from_bytes(encoded, 'little')

        if self.curve == "ed25519":
            y = integer & ((1 << 255) - 1)  # Extract lower 255 bits
            x_0 = (integer >> 255) & 1      # Extract MSB (LSB of x)
        else:  # ed448
            y = integer & ((1 << 447) - 1)
            x_0 = (integer >> 447) & 1

        # Ensure `y` is valid
        if y >= p:
            raise ValueError("Decoding failed: y >= p")

        # Compute `x^2` using curve equation
        y2 = (y * y) % p
        numerator = (y2 - 1) % p
        denominator = (d * y2 - a) % p

        # Compute modular inverse
        denom_inv = inverse_mod(denominator, p)

        # Compute `x^2 (mod p)`
        x2 = (numerator * denom_inv) % p

        # Ensure `x^2` is a quadratic residue
        if legendre_symbol(x2, p) != 1:
            raise ValueError(f"x^2 = {x2} is not a quadratic residue mod p")

        # Compute `x` using Tonelli-Shanks
        x = self.tonelli_shanks(x2, p)

        # Ensure `x` has the correct LSB
        if (x % 2) != x_0:
            x = p - x

        return x, y
    




    def key_generation(self):

        b = self.b

        # 1. Generate private key `d` (random `b`-bit string)
        d = os.urandom(b // 8)

        # 2. Compute hash `H(d)`
        hdigest = self.hash(d)

        # 3. Extract `hdigest1` (first half of `H(d)`)
        hdigest1 = bytearray(hdigest[:b // 8])  # First `b/8` bytes

        # 3.1 Modify `hdigest1` for Ed25519
        if self.curve == "ed25519":
            hdigest1[0] &= 0b11111000  # Set first 3 bits to 0
            hdigest1[-1] &= 0b01111111  # Set last bit to 0
            hdigest1[-1] |= 0b01000000  # Set second to last bit to 1

        # 3.2 Modify `hdigest1` for Ed448
        else:  # ed448
            hdigest1[0] &= 0b11111100  # Set first 2 bits to 0
            hdigest1[-1] = 0  # Set last octet to 0
            hdigest1[-2] |= 0b00000001  # Set last bit of second to last octet to 1

        # 4. Convert `hdigest1` to integer `s` (little-endian)
        self.s = int.from_bytes(hdigest1, byteorder="little")

        # 5. Compute public key `Q = [s]G`
        Q = self.G.mult(self.s)

        # Return (private key, public key)
        Q_encoded = self.encode_point(Q)

        return d, Q_encoded  
    

    def dom4(f, c):
        # "SigEd448" as a byte string in ASCII (8 bytes)
        sig_ed448 = b"SigEd448"

        # Convert the value f to a single byte (octet with value f)
        f_octet = bytes([f])

        # Convert the length of the string c to a 2-byte big-endian representation
        c_length = len(c).to_bytes(2, 'little')

        # Return the concatenation of "SigEd448" || f || len(c) || c
        return sig_ed448 + f_octet + c_length + c.encode('ascii')


    def sign(self, M, d, Q, context=b""):

        b = self.b  # Bit length (256 for Ed25519, 456 for Ed448)

        # Compute hash H(d)
        hdigest = self.hash(d)

        # Extract the second half of H(d)
        hdigest2 = hdigest[b // 8:]

        # Compute r
        if self.curve == "ed25519":
            r = hashlib.sha512(hdigest2 + M).digest()
            r = int.from_bytes(r, "little")
        else:  # ed448
            dom = self.dom4(0,context)
            r = hashlib.shake_256(dom + hdigest2 + M).digest(912//8)
            r = int.from_bytes(r, "little")

        # Compute point R = [r]G and encode it
        R = self.G.mult(r)
        R_encoded = self.encode_point(R)

        # Compute `S`
        if self.curve == "ed25519":
            h = hashlib.sha512(R_encoded + Q + M).digest()
            h = int.from_bytes(h, "little")
        else:  # ed448
            h = hashlib.shake_256(dom + R_encoded + Q + M).digest(912//8)
            h = int.from_bytes(h, "little")

        S = (r + h * self.s) % self.n
        S_encoded = int(S).to_bytes(self.b//8, 'little')

        # Return the signature R || S
        return R_encoded + S_encoded





    def verify(self, M, signature, Q, context=b""):

        # Step 1: Decode the signature and public key
        if len(signature) != 2 * (self.b // 8):
            return False  # Invalid signature length
        
        R_bytes = signature[:self.b // 8]
        S_bytes = signature[self.b // 8:]

        # Decode S as an integer
        s = int.from_bytes(S_bytes, byteorder='little')

        # Verify that S is in range [0, n)
        if s < 0 or s >= self.n:
            return False  # S out of range

        # Decode R as a point
        try:
            R_x, R_y = self.decode_point(R_bytes)
            R = ed(curve=self.E, x=R_x, y=R_y)
        except ValueError:
            return False  # R decoding failed

        # Decode Q as a point
        try:
            Q_x, Q_y = self.decode_point(Q)
            Q_point = ed(curve=self.E, x=Q_x, y=Q_y)
        except ValueError:
            return False  # Q decoding failed

        # Step 2: Form HashData = R || Q || M
        HashData = R_bytes + Q + M

        # Step 3: Compute digest
        if self.curve == "ed25519":
            # Step 3.1: For Ed25519
            digest = hashlib.sha512(HashData).digest()
        else:  # ed448
            # Step 3.2: For Ed448
            dom = b"SigEd448" + bytes([0, len(context)]) + context
            digest = hashlib.shake_256(dom + HashData).digest(912//8)

        # Interpret digest as little-endian integer t
        t = int.from_bytes(digest, byteorder='little')


        # Step 4: Check verification equation [2^c * S]G = [2^c]R + (2^c * t)Q
        # Calculate left side: [2^c * S]G
        left_side = self.G.mult(2 ** self.c * s)

        # Calculate right side: [2^c]R + (2^c * t)Q
        R_scaled = R.mult(2 ** self.c)
        Q_scaled = Q_point.mult(2 ** self.c * t)

        # Add the points
        right_side = R_scaled.copy()
        right_side.soma(Q_scaled)

        # Compare left and right sides
        return left_side.eq(right_side)




# Test class for EdDSA
def test_eddsa():
    # Initialize and test both curve types
    test_curves = ["ed25519"]
    
    for curve_name in test_curves:
        print(f"\n{'=' * 60}")
        print(f"Testing {curve_name}")
        print(f"{'=' * 60}")
        
        # Initialize the ECDSA object
        print("Initializing curve...")
        try:
            ecdsa = ECDSA(curve_name)
            print("✓ Successfully initialized curve")
        except Exception as e:
            print(f"✗ Failed to initialize curve: {e}")
            continue
        
        # Test key generation
        print("\nTesting key generation...")
        try:
            start_time = time.time()
            private_key, public_key = ecdsa.key_generation()
            end_time = time.time()
            
            print(f"✓ Key generation completed in {end_time - start_time:.4f} seconds")
            print(f"✓ Private key length: {len(private_key)} bytes")
            print(f"✓ Public key length: {len(public_key)} bytes")
        except Exception as e:
            print(f"✗ Key generation failed: {e}")
            continue
        
        # Test signing
        print("\nTesting message signing...")
        messages = [
            b"This is a test message",
            b"Another test with different content",
            b"A third test message with more data to sign"
        ]
        
        for i, message in enumerate(messages):
            try:
                start_time = time.time()
                signature = ecdsa.sign(message, private_key, public_key)
                end_time = time.time()
                
                print(f"✓ Message {i+1} signed in {end_time - start_time:.4f} seconds")
                print(f"✓ Signature length: {len(signature)} bytes")
                
                # Test verification
                start_time = time.time()
                verification_result = ecdsa.verify(message, signature, public_key)
                end_time = time.time()
                
                if verification_result:
                    print(f"✓ Signature {i+1} verified successfully in {end_time - start_time:.4f} seconds")
                else:
                    print(f"✗ Signature {i+1} verification failed")
            except Exception as e:
                print(f"✗ Signing or verification failed for message {i+1}: {e}")
        
        # Test with context (especially for Ed448)
        if curve_name == "ed448":
            print("\nTesting with context (Ed448)...")
            context = b"test-context"
            try:
                message = b"Message with context"
                signature = ecdsa.sign(message, private_key, public_key, context)
                verification_result = ecdsa.verify(message, signature, public_key, context)
                
                if verification_result:
                    print(f"✓ Signature with context verified successfully")
                else:
                    print(f"✗ Signature with context verification failed")
            except Exception as e:
                print(f"✗ Signing or verification with context failed: {e}")
        
        # Test with invalid signature
        print("\nTesting with invalid signature...")
        try:
            message = b"Original message"
            signature = ecdsa.sign(message, private_key, public_key)
            
            # Modify the signature to make it invalid
            modified_signature = bytearray(signature)
            modified_signature[0] ^= 0xFF  # Flip bits in the first byte
            modified_signature = bytes(modified_signature)
            
            verification_result = ecdsa.verify(message, modified_signature, public_key)
            
            if not verification_result:
                print(f"✓ Invalid signature correctly rejected")
            else:
                print(f"✗ Invalid signature incorrectly accepted")
        except Exception as e:
            print(f"✗ Invalid signature test failed: {e}")
        
        # Test with wrong message
        print("\nTesting with wrong message...")
        try:
            original_message = b"Original message"
            signature = ecdsa.sign(original_message, private_key, public_key)
            
            wrong_message = b"Wrong message"
            verification_result = ecdsa.verify(wrong_message, signature, public_key)
            
            if not verification_result:
                print(f"✓ Signature for wrong message correctly rejected")
            else:
                print(f"✗ Signature for wrong message incorrectly accepted")
        except Exception as e:
            print(f"✗ Wrong message test failed: {e}")


# Run the test if this script is executed directly
if __name__ == "__main__":
    print("Starting EdDSA implementation test...")
    test_eddsa()