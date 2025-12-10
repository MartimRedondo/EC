# EC-TP4 – Implementação do Algoritmo Dilithium em SageMath

## 1. O que está implementado com sucesso

### Geração de Chaves (keygen)

- Gera a chave pública `pk` e privada `sk`, com:
    - `rho`: semente para gerar a matriz pública `A`;
    - `s1`, `s2`: vetores secretos;
    - `t = A*s1 + s2` (em NTT);
    - `t0` e `t1`: divisão em parte baixa/alta;
    - Empacotamento final de `pk` e `sk`.

### Assinatura (sign)

- Implementada com **rejection sampling**:
    - Gera `y`, calcula `w = A*y`, extrai `w1 = HighBits(w)`;
    - Gera `c = H(mu || w1)`, calcula `z = y + c·s1`;
    - Verifica o limite de norma;
    - Se válido, empacota assinatura `sig = c || z || h` (onde `h` ainda é um placeholder);
    - Guarda dados de debug (`w1`, `mu`, `A`, etc.).

### Verificação (verify)

- Desempacota `z`, `c`, `pk`;
- Recalcula `mu = H(H(pk) || mensagem)`;
- Reconstrói `w1' = HighBits(A*z - c*t1)`;
- Compara `c' = H(mu || w1')` com `c` da assinatura.

### Testes Unitários

- `test_dilithium_sage.py` cobre:
    - Assinaturas válidas/inválidas;
    - Mensagens alteradas e assinaturas corrompidas;
    - Verificação de aleatoriedade e comprimento.

- `test_polynomials_sage.py` e `test_modules_sage.py`:
    - Testam aritmética de polinómios e matrizes;
    - NTT, packing, amostragem, etc.
    - Todos os testes passam com sucesso.

---

## 2. Problemas Encontrados e Funcionalidades Incompletas

### Verificação com erro

- A função `verify()` rejeita a assinatura por divergência entre os `w1_bytes` recomputados e os originais:
    - `A`, `mu` e `c` coincidem entre `sign()` e `verify()`;
    - Mas `w1` reconstruído é diferente;
    - Isso leva a `c' != c` → assinatura inválida.

#### Causas prováveis:
- Erro no `_reconstruct_w1()` (ex: falta de `from_ntt()` antes do `power_2_round`);
- Arredondamento ou mod q aplicado incorretamente;
- Erro no packing dos `t1` (coeficientes fora de intervalo).

### `h` (Hints) não implementados

- Na assinatura:
    - `h_bytes = b"\x00" * self.omega` é apenas um placeholder;
    - `make_hint()` e `use_hint()` não são usados;
    - Isto afeta a segurança e a compacidade da assinatura.

### Verificação de limites de `z` incompleta

- Apenas verificado na assinatura (sign), não na verificação (verify);
- Normalmente deve-se testar se `||z||∞ < γ1 - β` também no `verify()`.

---

## 3. Adaptação para SageMath – Módulos Base

### polynomials_sage.py

- Define a estrutura `PolynomialDilithium` e `PolynomialRingDilithium`;
- Suporte para:
    - Aritmética polinomial;
    - NTT (e inversa);
    - Empacotamento/Desempacotamento (`bit_pack_t0`, `bit_unpack_t0`);
    - `power_2_round`, `decompose`, `make_hint`, `use_hint`;
    - Amostragem (`sample_in_ball`, `rejection_sample_ntt_poly`, etc.);
- Todos os testes passam com sucesso.

### modules_sage.py

- Define estruturas de matrizes (`MatrixDilithium`, `ModuleDilithium`);
- Operações:
    - Soma, multiplicação, transposição, conversão NTT;
    - Geração de matriz `A` com seed e vetores NTT;
    - Desempacotamento de polinómios tipo `t1`, `t0`, `z`, etc.
- Testes locais passaram com sucesso.

---

## 4. Explicação do ficheiro dilithium_sage.py

### Estrutura Geral

- Classe principal: `Dilithium`;
- Métodos:
    - `keygen()`: gera chaves;
    - `sign(m, sk, pk)`: assina a mensagem `m`;
    - `verify(m, sig, pk)`: verifica a assinatura;
    - `_generate_challenge`, `_reconstruct_w1`, packing/unpacking, etc.

### Processo Resumido:

- **keygen()**:
    - Gera `A`, `s1`, `s2`;
    - Calcula `t = A*s1 + s2`;
    - Divide em `t0`, `t1` e empacota as chaves.

- **sign()**:
    - Gera digest `mu`;
    - Executa rejection sampling para encontrar `z`;
    - Assinatura = `c || z || h`.

- **verify()**:
    - Desempacota `c`, `z`, `pk`;
    - Recalcula `mu`, reconstrói `w1'`;
    - Compara `c'` com `c`.

---

## 5. Conclusão e Trabalhos Futuros

Apesar da base estar sólida e os testes dos módulos passarem com sucesso, a assinatura ainda falha na verificação por divergência nos `w1_bytes`. Acreditamos que o erro está isolado na função `_reconstruct_w1()` e poderá ser resolvido com:

1. **Revisão fina da transformação NTT e `power_2_round`**;
2. **Comparação linha-a-linha de `w1` (versão `sign()` vs `verify()`)**;
3. **Integração das funções `make_hint` e `use_hint` para completar a verificação leve.**

Estamos disponíveis para apresentar esta solução se for necessário para avaliação.

