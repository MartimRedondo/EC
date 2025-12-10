"""
Versão instrumentada do código de decapsulamento BIKE com debug detalhado
para identificar a causa da não coincidência dos segredos compartilhados.
"""

import hashlib
import os
import time
import logging
import numpy as np
from typing import Tuple, List, Optional, Dict, Any
from sage.all import *

# Configuração de logging detalhado
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('bike_debug.log')
    ]
)
logger = logging.getLogger("BIKE_DEBUG")

# Função para comparar valores e registrar diferenças
def compare_values(name, value1, value2, context=""):
    """Compara dois valores e registra diferenças detalhadas"""
    if isinstance(value1, list) and isinstance(value2, list):
        if len(value1) != len(value2):
            logger.error(f"{context} - {name}: Comprimentos diferentes: {len(value1)} vs {len(value2)}")
            return False
        
        differences = sum(1 for i in range(len(value1)) if value1[i] != value2[i])
        if differences > 0:
            logger.error(f"{context} - {name}: {differences}/{len(value1)} elementos diferentes")
            # Mostrar as primeiras diferenças
            for i in range(min(len(value1), len(value2))):
                if value1[i] != value2[i]:
                    logger.error(f"  Posição {i}: {value1[i]} vs {value2[i]}")
                    if i >= 10:  # Limitar a 10 diferenças
                        logger.error("  ... mais diferenças omitidas")
                        break
            return False
        return True
    else:
        is_equal = value1 == value2
        if not is_equal:
            logger.error(f"{context} - {name}: Valores diferentes: {value1} vs {value2}")
        return is_equal

# Função para registrar estado de polinômios
def log_polynomial_state(name, poly, context=""):
    """Registra informações detalhadas sobre um polinômio"""
    if poly is None:
        logger.debug(f"{context} - {name}: None")
        return
    
    coeffs = poly.list()
    weight = sum(1 for c in coeffs if c != 0)
    degree = poly.degree() if poly != 0 else -1
    
    logger.debug(f"{context} - {name}: grau={degree}, peso={weight}, primeiros coefs={coeffs[:10]}")
    
    # Registrar posições dos coeficientes não-zero (útil para polinômios esparsos)
    if weight > 0 and weight <= 200:  # Limitar para polinômios esparsos
        non_zero_pos = [i for i, c in enumerate(coeffs) if c != 0]
        logger.debug(f"{context} - {name} posições não-zero: {non_zero_pos}")

# Função para registrar estado de vetores
def log_vector_state(name, vec, context=""):
    """Registra informações detalhadas sobre um vetor"""
    if vec is None:
        logger.debug(f"{context} - {name}: None")
        return
    
    weight = sum(1 for bit in vec if bit != 0)
    
    logger.debug(f"{context} - {name}: tamanho={len(vec)}, peso={weight}, primeiros bits={list(vec[:20])}")
    
    # Registrar posições dos bits 1 (útil para vetores esparsos)
    if weight > 0 and weight <= 200:  # Limitar para vetores esparsos
        ones_pos = [i for i, bit in enumerate(vec) if bit != 0]
        logger.debug(f"{context} - {name} posições de 1s: {ones_pos}")

# Função para instrumentar o decapsulamento
def instrumented_decaps(ct, sk,pk, original_decaps_func, params, hash_functions, bgf_decoder):
    """
    Versão instrumentada da função de decapsulamento com debug detalhado
    para identificar onde ocorre a divergência.
    """
    logger.info("=" * 80)
    logger.info("INICIANDO DECAPSULAMENTO INSTRUMENTADO")
    logger.info("=" * 80)
    
    # Registrar parâmetros
    logger.debug(f"Parâmetros: R_BITS={params.R_BITS}, DV={params.DV}, T1={params.T1}")
    
    # Registrar estado das chaves
    log_polynomial_state("sk.h0", sk.h0, "CHAVES")
    log_polynomial_state("sk.h1", sk.h1, "CHAVES")
    log_vector_state("sk.sigma", sk.sigma, "CHAVES")
    
    # Registrar estado do ciphertext
    log_polynomial_state("ct.c0", ct.c0, "CIPHERTEXT")
    log_vector_state("ct.c1", ct.c1, "CIPHERTEXT")
    
    try:
        # ETAPA 1: Cálculo da síndrome
        logger.info("ETAPA 1: Calculando síndrome")
        start_time = time.time()
        
        # Registrar produto antes da redução
        prod_raw = ct.c0 * sk.h0
        log_polynomial_state("prod_raw (c0*h0)", prod_raw, "SÍNDROME")
        
        # Calcular síndrome (produto módulo polinômio)
        mod_poly_r = params.mod_poly_r if hasattr(params, 'mod_poly_r') else (params.x**params.R_BITS - 1)
        prod = (ct.c0 * sk.h0) % mod_poly_r
        
        # Extrair coeficientes e converter para vetor
        coeffs = prod.list()
        if len(coeffs) < params.R_BITS:
            coeffs.extend([0] * (params.R_BITS - len(coeffs)))
        
        syndrome = params.vector_class([params.field_class(int(c)) for c in coeffs[:params.R_BITS]])
        
        log_vector_state("syndrome", syndrome, "SÍNDROME")
        logger.debug(f"Tempo de cálculo da síndrome: {time.time() - start_time:.6f}s")
        
        # ETAPA 2: Decodificação BGF
        logger.info("ETAPA 2: Decodificando com BGF")
        start_time = time.time()
        
        # Chamar o decodificador BGF
        e_prime_poly = bgf_decoder.decode(syndrome, sk.h0, sk.h1)
        
        log_polynomial_state("e_prime_poly", e_prime_poly, "BGF")
        logger.debug(f"Tempo de decodificação BGF: {time.time() - start_time:.6f}s")
        
        # ETAPA 3: Recuperação da mensagem
        logger.info("ETAPA 3: Recuperando mensagem")
        start_time = time.time()
        
        # Calcular L'
        L_prime = hash_functions.function_L(e_prime_poly)
        log_vector_state("L_prime", L_prime, "RECUPERAÇÃO")
        
        # Recuperar m' = c1 ⊕ L'
        m_prime = params.xor_func(ct.c1, L_prime)
        log_vector_state("m_prime", m_prime, "RECUPERAÇÃO")
        
        logger.debug(f"Tempo de recuperação da mensagem: {time.time() - start_time:.6f}s")
        
        # ETAPA 4: Verificação de consistência
        logger.info("ETAPA 4: Verificando consistência")
        start_time = time.time()
        
        # Recomputar vetor de erro a partir de m'
        e_check = hash_functions.function_H_deterministic(m_prime)
        log_polynomial_state("e_check", e_check, "VERIFICAÇÃO")
        
        # Dividir o vetor de erro
        e0_check, e1_check = params.split_func(e_check)
        log_polynomial_state("e0_check", e0_check, "VERIFICAÇÃO")
        log_polynomial_state("e1_check", e1_check, "VERIFICAÇÃO")
        
        # Reconstruir h
        h0_inv = params.inverse_func(sk.h0)
        log_polynomial_state("h0_inv", h0_inv, "VERIFICAÇÃO")
        
        h_reconstructed = (sk.h1 * h0_inv) % mod_poly_r
        log_polynomial_state("h_reconstructed", h_reconstructed, "VERIFICAÇÃO")
        
        # Recomputar c0
        temp_check = (e1_check * pk.h) % params.mod_poly_r
        log_polynomial_state("temp_check", temp_check, "VERIFICAÇÃO")
        
        c0_check = (e0_check + temp_check) % params.mod_poly_r
        log_polynomial_state("c0_check", c0_check, "VERIFICAÇÃO")
        
        # Verificar consistência
        diff = (ct.c0 - c0_check) % mod_poly_r
        is_consistent = (diff == 0)
        
        logger.info(f"Consistência: {is_consistent}")
        logger.debug(f"Tempo de verificação de consistência: {time.time() - start_time:.6f}s")
        
        # ETAPA 5: Derivação da chave
        logger.info("ETAPA 5: Derivando chave compartilhada")
        start_time = time.time()
        
        # Derivar chave baseada na consistência
        if is_consistent:
            logger.info("Usando mensagem recuperada para derivar chave")
            k = hash_functions.function_K(m_prime, ct.c0, ct.c1)
        else:
            logger.warning("Inconsistência detectada - usando fallback sigma")
            k = hash_functions.function_K(sk.sigma, ct.c0, ct.c1)
        
        log_vector_state("k (chave derivada)", k, "CHAVE")
        logger.debug(f"Tempo de derivação da chave: {time.time() - start_time:.6f}s")
        
        # Criar objeto de segredo compartilhado
        ss = params.shared_secret_class(k)
        
        logger.info("=" * 80)
        logger.info("DECAPSULAMENTO INSTRUMENTADO CONCLUÍDO")
        logger.info("=" * 80)
        
        return ss
        
    except Exception as e:
        logger.exception(f"ERRO NO DECAPSULAMENTO INSTRUMENTADO: {e}")
        raise

# Função para instrumentar o encapsulamento
def instrumented_encaps(pk, original_encaps_func, params, hash_functions):
    """
    Versão instrumentada da função de encapsulamento com debug detalhado
    para comparação com o decapsulamento.
    """
    logger.info("=" * 80)
    logger.info("INICIANDO ENCAPSULAMENTO INSTRUMENTADO")
    logger.info("=" * 80)
    
    # Registrar parâmetros
    logger.debug(f"Parâmetros: R_BITS={params.R_BITS}, DV={params.DV}, T1={params.T1}")
    
    # Registrar estado da chave pública
    log_polynomial_state("pk.h", pk.h, "CHAVE_PÚBLICA")
    
    try:
        # ETAPA 1: Geração da mensagem
        logger.info("ETAPA 1: Gerando mensagem aleatória")
        start_time = time.time()
        
        # Gerar mensagem aleatória
        m = params.secure_random_func(params.ELL_BITS)
        log_vector_state("m", m, "MENSAGEM")
        
        logger.debug(f"Tempo de geração da mensagem: {time.time() - start_time:.6f}s")
        
        # ETAPA 2: Geração do vetor de erro
        logger.info("ETAPA 2: Gerando vetor de erro")
        start_time = time.time()
        
        # Gerar vetor de erro
        e_poly = hash_functions.function_H_deterministic(m)
        log_polynomial_state("e_poly", e_poly, "ERRO")
        
        # Dividir o vetor de erro
        e0_poly, e1_poly = params.split_func(e_poly)
        log_polynomial_state("e0_poly", e0_poly, "ERRO")
        log_polynomial_state("e1_poly", e1_poly, "ERRO")
        
        logger.debug(f"Tempo de geração do vetor de erro: {time.time() - start_time:.6f}s")
        
        # ETAPA 3: Cálculo do ciphertext
        logger.info("ETAPA 3: Calculando ciphertext")
        start_time = time.time()
        
        # Calcular c0
        mod_poly_r = params.mod_poly_r if hasattr(params, 'mod_poly_r') else (params.x**params.R_BITS - 1)
        temp = (e1_poly * pk.h) % mod_poly_r
        log_polynomial_state("temp (e1*h)", temp, "CIPHERTEXT")
        
        c0 = (e0_poly + temp) % mod_poly_r
        log_polynomial_state("c0", c0, "CIPHERTEXT")
        
        # Calcular c1
        L_vec = hash_functions.function_L(e_poly)
        log_vector_state("L_vec", L_vec, "CIPHERTEXT")
        
        c1 = params.xor_func(L_vec, m)
        log_vector_state("c1", c1, "CIPHERTEXT")
        
        logger.debug(f"Tempo de cálculo do ciphertext: {time.time() - start_time:.6f}s")
        
        # ETAPA 4: Derivação da chave
        logger.info("ETAPA 4: Derivando chave compartilhada")
        start_time = time.time()
        
        # Derivar segredo compartilhado
        ss = hash_functions.function_K(m, c0, c1)
        log_vector_state("ss (chave derivada)", ss, "CHAVE")
        
        logger.debug(f"Tempo de derivação da chave: {time.time() - start_time:.6f}s")
        
        # Criar objetos de retorno
        ct = params.ciphertext_class(c0, c1)
        shared_secret = params.shared_secret_class(ss)
        
        logger.info("=" * 80)
        logger.info("ENCAPSULAMENTO INSTRUMENTADO CONCLUÍDO")
        logger.info("=" * 80)
        
        return ct, shared_secret
        
    except Exception as e:
        logger.exception(f"ERRO NO ENCAPSULAMENTO INSTRUMENTADO: {e}")
        raise

# Função para comparar segredos compartilhados
def compare_shared_secrets(ss_enc, ss_dec):
    """Compara detalhadamente os segredos compartilhados"""
    logger.info("=" * 80)
    logger.info("COMPARAÇÃO DE SEGREDOS COMPARTILHADOS")
    logger.info("=" * 80)
    
    enc_raw = list(ss_enc.raw)
    dec_raw = list(ss_dec.raw)
    
    match = compare_values("Segredos compartilhados", enc_raw, dec_raw, "COMPARAÇÃO")
    
    if match:
        logger.info("✅ SUCESSO: Segredos compartilhados coincidem!")
    else:
        logger.error("❌ FALHA: Segredos compartilhados não coincidem")
        
        # Análise adicional
        if len(enc_raw) == len(dec_raw):
            diff_count = sum(1 for i in range(len(enc_raw)) if enc_raw[i] != dec_raw[i])
            diff_percent = (diff_count / len(enc_raw)) * 100
            logger.error(f"Diferenças: {diff_count}/{len(enc_raw)} bits ({diff_percent:.2f}%)")
            
            # Mostrar as primeiras diferenças
            for i in range(min(20, len(enc_raw))):
                if enc_raw[i] != dec_raw[i]:
                    logger.error(f"  Bit {i}: enc={enc_raw[i]}, dec={dec_raw[i]}")
    
    return match

def run_instrumented_test():
    """Executa teste com instrumentação detalhada"""
    
    # Importar as funções e classes necessárias
    from seu_notebook import (
        generate_keypair, encaps, decaps_improved,
        HashFunctions, improved_bgf_decoder,
        PARAMS, F, R, x, mod_poly_r, mod_poly_n,
        split_vector_improved, xor_vectors_improved,
        safe_polynomial_inverse, vector,
        SecretKey, PublicKey, Ciphertext, SharedSecret
    )
    
    # Criar objeto de parâmetros instrumentados
    params = type('Params', (), {
        'R_BITS': PARAMS.R_BITS,
        'N_BITS': PARAMS.N_BITS,
        'DV': PARAMS.DV,
        'T1': PARAMS.T1,
        'ELL_BITS': PARAMS.ELL_BITS,
        'mod_poly_r': mod_poly_r,
        'mod_poly_n': mod_poly_n,
        'x': x,
        'F': F,
        'R': R,
        'field_class': F,
        'vector_class': vector,
        'ciphertext_class': Ciphertext,
        'shared_secret_class': SharedSecret,
        'secure_random_func': lambda bits: SecureRandom.secure_vector(bits),
        'split_func': split_vector_improved,
        'xor_func': xor_vectors_improved,
        'inverse_func': safe_polynomial_inverse
    })()
    
    try:
        # Gerar chaves
        logger.info("Gerando par de chaves...")
        pk, sk = generate_keypair()
        
        # Log detalhado das chaves
        logger.debug(f"Peso de h0: {sum(1 for c in sk.h0.list() if c != 0)}")
        logger.debug(f"Peso de h1: {sum(1 for c in sk.h1.list() if c != 0)}")
        
        # Encapsulamento instrumentado
        ct, ss_enc = instrumented_encaps(pk, encaps, params, HashFunctions)
        
        # Decapsulamento instrumentado  
        ss_dec = instrumented_decaps(ct, sk, decaps_improved, params, 
                                     HashFunctions, improved_bgf_decoder, pk)
        
        # Comparar resultados
        match = compare_shared_secrets(ss_enc, ss_dec)
        
        return match
        
    except Exception as e:
        logger.exception(f"Erro no teste instrumentado: {e}")
        raise