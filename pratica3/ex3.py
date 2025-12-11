from sage.all import *
import hashlib
import secrets
from enum import Enum
import random

# ======= Parte 1: Implementação do Circuito Booleano =======
class GateType(Enum):
    AND = 0
    XOR = 1

class Gate:
    def __init__(self, gate_type, input1, input2, output):
        self.gate_type = gate_type  # AND ou XOR
        self.input1 = input1        # Índice da primeira entrada
        self.input2 = input2        # Índice da segunda entrada
        self.output = output        # Índice da saída
        
    def __str__(self):
        gate_name = "AND" if self.gate_type == GateType.AND else "XOR"
        return f"{gate_name}: ({self.input1}, {self.input2}) -> {self.output}"

class Circuit:
    def __init__(self, n_inputs):
        self.n_inputs = n_inputs      # Número de entradas
        self.gates = []               # Lista de portas
        self.wire_count = n_inputs    # Contador de fios (começa com os inputs)
        
    def add_gate(self, gate_type, input1, input2):
        """Adiciona uma nova porta ao circuito"""
        new_gate = Gate(gate_type, input1, input2, self.wire_count)
        self.gates.append(new_gate)
        self.wire_count += 1
        return self.wire_count - 1    # Retorna o índice do fio de saída
    
    def get_output_wire(self):
        """Retorna o índice do fio de saída (último adicionado)"""
        if not self.gates:
            return 0  # Se não há portas, retorna a primeira entrada
        return self.wire_count - 1
    
    def __str__(self):
        result = f"Circuit com {self.n_inputs} entradas e {len(self.gates)} portas\n"
        for i, gate in enumerate(self.gates):
            result += f"Porta {i}: {gate}\n"
        result += f"Fio de saída: {self.get_output_wire()}"
        return result

# ======= Parte 2: Função de expansão de saída (XOF) usando SageMath =======
def create_xof(seed, lambda_param=256):
    """
    Cria um XOF a partir de uma seed usando SageMath
    
    Args:
        seed: Semente para o gerador (inteiro ou bytes)
        lambda_param: Tamanho da segurança em bits
    
    Returns:
        Uma função que pode ser usada para extrair bits do XOF
    """
    # Convertemos a seed para bytes se for um inteiro
    if isinstance(seed, Integer) or isinstance(seed, int):
        # Em SageMath, usamos Integer para trabalhar com precisão arbitrária
        seed_bytes = int(seed).to_bytes((int(seed).bit_length() + 7) // 8, byteorder='big')
    else:
        seed_bytes = seed
        
    # Inicializa o SHAKE256 com a seed
    shake = hashlib.shake_256(seed_bytes)
    
    # Função para extrair bits do XOF
    def get_bits(n_bits):
        # Obter bytes suficientes para n_bits
        n_bytes = (n_bits + 7) // 8
        bytes_data = shake.digest(n_bytes)
        
        # Converter bytes para um inteiro usando SageMath
        value = Integer(int.from_bytes(bytes_data, byteorder='big'))
        
        # Máscara para obter apenas n_bits
        mask = (1 << n_bits) - 1
        return value & mask
    
    return get_bits

# ======= Parte 3: Geração do Circuito Booleano =======
def generate_polynomial_circuit(seed, n, polynomial_degree=2):
    """
    Gera um circuito booleano de tamanho poly(n) a partir de uma seed
    
    Args:
        seed: Semente para o gerador de números aleatórios
        n: Número de entradas do circuito
        polynomial_degree: Grau do polinômio que determina o tamanho do circuito
    
    Returns:
        Um circuito booleano
    """
    # Calculamos o número de portas usando um polinômio
    # Usando SageMath para calcular n^polynomial_degree
    n_gates = Integer(n)^polynomial_degree + 2*n
    
    # Cria o XOF
    xof = create_xof(seed)    

    # Cria o circuito
    circuit = Circuit(n)
    
    # Adiciona portas ao circuito
    available_wires = list(range(n))  # Inicialmente, apenas as entradas estão disponíveis
    
    for _ in range(n_gates):
        # Escolhe o tipo da porta (AND ou XOR)
        gate_type = GateType.AND if xof(1) == 0 else GateType.XOR
        
        # Escolhe as entradas da porta entre os fios disponíveis
        input1_idx = xof(max(8, len(available_wires).bit_length() + 1))
        input1 = available_wires[input1_idx % len(available_wires)]
        
        input2_idx = xof(max(8, len(available_wires).bit_length() + 1))
        input2 = available_wires[input2_idx % len(available_wires)]
        
        # Adiciona a porta ao circuito
        output_wire = circuit.add_gate(gate_type, input1, input2)
        
        # Adiciona o fio de saída à lista de fios disponíveis
        available_wires.append(output_wire)
    
    # Garantir que o último fio é a saída do circuito
    # Se necessário, adicionar uma porta extra ligando o último fio de saída a algo
    if len(circuit.gates) > 0:
        last_output = circuit.get_output_wire()
        
        # Se ainda temos entradas disponíveis, criar uma porta adicional
        input1_idx = xof(max(8, len(available_wires).bit_length() + 1))
        input1 = available_wires[input1_idx % len(available_wires)]
        
        # A última porta será do tipo escolhido aleatoriamente 
        final_gate_type = GateType.AND if xof(1) == 0 else GateType.XOR
        circuit.add_gate(final_gate_type, last_output, input1)
    
    return circuit

def evaluate_circuit(circuit, inputs):
    """
    Avalia um circuito booleano com os inputs fornecidos
    
    Args:
        circuit: O circuito a ser avaliado
        inputs: Uma lista de booleanos com os valores de entrada
    
    Returns:
        O valor de saída do circuito (0 ou 1)
    """
    if len(inputs) != circuit.n_inputs:
        raise ValueError(f"Expected {circuit.n_inputs} inputs, got {len(inputs)}")
    
    # Inicializa os valores dos fios com os inputs
    wire_values = list(inputs)
    
    # Avalia cada porta em ordem topológica
    for gate in circuit.gates:
        input1_value = wire_values[gate.input1]
        input2_value = wire_values[gate.input2]
        
        # Calcula o valor de saída com base no tipo da porta
        # Usando operações em GF(2) do SageMath
        if gate.gate_type == GateType.AND:
            output_value = GF(2)(input1_value) * GF(2)(input2_value)
        else:  # XOR
            output_value = GF(2)(input1_value) + GF(2)(input2_value)
        
        # Armazena o valor no fio de saída
        wire_values.append(int(output_value))
    
    # Retorna o valor do último fio (saída do circuito)
    return wire_values[-1]

# ======= Parte 4: Implementação do Protocolo MPC in the Head =======
class ViewShare:
    def __init__(self, w_share, z_share, m_share):
        self.w_share = w_share  # Shares dos inputs
        self.z_share = z_share  # Bits de correlated randomness
        self.m_share = m_share  # Mensagens

    def __str__(self):
        return f"ViewShare(w_share={self.w_share}, z_share={self.z_share}, m_share={self.m_share})"

class MPCitHProver:
    def __init__(self, circuit, security_param=256):
        self.circuit = circuit
        self.security_param = security_param
        self.F = GF(2)  # Campo finito binário do SageMath
    
    def share_secret(self, secret, num_shares=3):
        """
        Divide um segredo em shares usando um esquema 2-out-of-3
        
        Args:
            secret: O segredo a ser compartilhado (um bit)
            num_shares: Número de shares a serem gerados
        
        Returns:
            Lista de shares
        """
        # Usar SageMath para trabalhar em GF(2)
        secret_gf2 = self.F(secret)
        shares = []
        
        # Gerar aleatoriamente n-1 shares
        for _ in range(num_shares - 1):
            shares.append(int(self.F.random_element()))
        
        # Calcular o último share para que a soma seja igual ao segredo
        last_share = secret_gf2
        for share in shares:
            last_share += self.F(share)
        
        shares.append(int(last_share))
        return shares
    
    def generate_correlated_randomness(self, seed, sid):
        """
        Gera 'correlated randomness' para os participantes
        
        Args:
            seed: Seed para o gerador
            sid: Identificador da sessão
        
        Returns:
            Lista de CRs para cada participante
        """
        # Criar XOF a partir da seed e sid
        combined_seed = str(seed) + str(sid)
        xof = create_xof(combined_seed.encode())
        
        # Gerar CRs para cada participante
        num_and_gates = sum(1 for gate in self.circuit.gates if gate.gate_type == GateType.AND)
        cr_bits = []
        
        # Gerar bits aleatórios para os primeiros dois participantes
        for i in range(2):
            participant_bits = [xof(1) for _ in range(num_and_gates)]
            cr_bits.append(participant_bits)
        
        # Calcular bits para o terceiro participante
        third_participant_bits = []
        for i in range(num_and_gates):
            # XOR dos bits correspondentes para satisfazer ζ0 ⊕ ζ1 ⊕ ζ2 = 0
            bit = (cr_bits[0][i] + cr_bits[1][i]) % 2
            third_participant_bits.append(bit)
        
        cr_bits.append(third_participant_bits)
        return cr_bits
    
    def prove(self, witness, sid):
        """
        Gera prova para o circuito
        
        Args:
            witness: Valores das entradas que satisfazem o circuito
            sid: Identificador da sessão
        
        Returns:
            Views dos três participantes
        """
        if len(witness) != self.circuit.n_inputs:
            raise ValueError(f"Expected {self.circuit.n_inputs} inputs, got {len(witness)}")
        
        # Verificar se o witness satisfaz o circuito
        if evaluate_circuit(self.circuit, witness) != 1:
            raise ValueError("Witness não satisfaz o circuito")
        
        # Gerar seeds aleatórias para cada participante
        seeds = [secrets.randbits(self.security_param) for _ in range(3)]
        
        # Compartilhar inputs entre os participantes
        input_shares = []
        for input_bit in witness:
            shares = self.share_secret(input_bit)
            input_shares.append(shares)
        
        # Transpor para obter as shares de cada participante
        participant_input_shares = list(map(list, zip(*input_shares)))
        
        # Gerar correlated randomness para as portas AND
        cr_bits = self.generate_correlated_randomness(sum(seeds), sid)
        
        # Simular a execução do protocolo MPC
        views = self.simulate_mpc(participant_input_shares, cr_bits)
        
        return views
    
    def simulate_mpc(self, participant_input_shares, cr_bits):
        """
        Simula a execução do protocolo MPC
        
        Args:
            participant_input_shares: Shares dos inputs para cada participante
            cr_bits: Bits de correlated randomness para cada participante
        
        Returns:
            Views dos três participantes
        """
        # Inicializar shares dos wires para cada participante
        wire_shares = [participant_input_shares[i][:] for i in range(3)]
        
        # Inicializar mensagens para gates AND
        messages = [[] for _ in range(3)]
        
        # Processar cada porta em ordem topológica
        and_gate_idx = 0
        for gate in self.circuit.gates:
            input1_idx = gate.input1
            input2_idx = gate.input2
            
            for p_idx in range(3):
                # Para cada participante, calcular sua share para esta porta
                input1_share = wire_shares[p_idx][input1_idx]
                input2_share = wire_shares[p_idx][input2_idx]
                
                if gate.gate_type == GateType.XOR:
                    # Para XOR, simplesmente somamos as shares em GF(2)
                    output_share = (input1_share + input2_share) % 2
                else:  # AND gate
                    # Para AND, usamos o protocolo específico
                    zeta = cr_bits[p_idx][and_gate_idx]
                    
                    # Calcular r_i
                    prev_idx = (p_idx - 1) % 3
                    prev_zeta = cr_bits[prev_idx][and_gate_idx]
                    
                    # O cálculo de r_i envolve multiplicações das shares 
                    # e o bit de correlated randomness
                    r_i = (input1_share * input2_share) ^ zeta
                    
                    # Armazenar r_i como mensagem enviada ao próximo participante
                    next_idx = (p_idx + 1) % 3
                    messages[p_idx].append(r_i)
                    
                    # O output share para AND é calculado usando r_(i-1) e r_i
                    if p_idx == 0:
                        # Para o primeiro participante, ainda não temos r_(i-1)
                        # mas na simulação podemos calculá-lo
                        r_prev = messages[prev_idx][-1] if messages[prev_idx] else 0
                    else:
                        r_prev = messages[prev_idx][-1]
                    
                    output_share = r_prev ^ r_i
                    
                    # Incrementar o índice da porta AND
                    if p_idx == 2:  # Apenas incrementar uma vez por gate
                        and_gate_idx += 1
                
                # Armazenar o output share
                wire_shares[p_idx].append(output_share)
        
        # Construir views para cada participante
        views = []
        for p_idx in range(3):
            w_share = participant_input_shares[p_idx]
            z_share = cr_bits[p_idx]
            m_share = messages[p_idx]
            views.append(ViewShare(w_share, z_share, m_share))
        
        return views

class ObliviousTransfer:
    def __init__(self, security_param=256):
        self.security_param = security_param
        
    def setup(self):
        """Configuração inicial do protocolo OT"""
        # Definir intervalo seguro com base no parâmetro de segurança
        min_p = 2 ** (self.security_param - 1)
        max_p = 2 ** self.security_param - 1

        # Garante que o intervalo tem pelo menos alguns valores válidos
        if max_p <= min_p:
            max_p = min_p + 1000  # ampliar artificialmente para garantir existência de primos

        # Tentar gerar um primo neste intervalo
        p = random_prime(max_p, lbound=min_p)
        self.field = GF(p)

        # Gerar um gerador do grupo
        self.g = self.field.multiplicative_generator()

        return p, self.g
    
    def sender_initial(self, messages):
        """
        Fase inicial do sender no protocolo OT 2-out-of-3
        
        Args:
            messages: Lista de 3 mensagens a serem transmitidas
        
        Returns:
            Dados públicos do sender
        """
        # Gerar chaves aleatórias
        self.sender_keys = [self.field.random_element() for _ in range(3)]
        
        # Calcular chaves públicas
        sender_public = [power_mod(int(self.g), int(key), int(self.field.order())) for key in self.sender_keys]
        
        # Armazenar mensagens para uso posterior
        self.messages = messages
        
        return sender_public

    def receiver_choose(self, sender_public, choices):
        """
        Receiver escolhe duas das três mensagens

        Args:
            sender_public: Chaves públicas do sender
            choices: Índices das duas mensagens escolhidas

        Returns:
            Dados enviados ao sender
        """
        if len(choices) != 2 or choices[0] == choices[1]:
            raise ValueError("Deve escolher exatamente dois índices distintos")
        
        # Armazenar escolhas
        self.choices = sorted(choices)
        
        # Gerar valores aleatórios para cada escolha
        self.receiver_values = [self.field.random_element() for _ in range(2)]
        
        # Calcular valores a serem enviados ao sender
        receiver_public = []
        for i, choice in enumerate(self.choices):
            # Obter o módulo correto diretamente do campo
            order = self.field.order()
            
            # Certificar que o 'order' seja um número inteiro
            n = int(order)  # Convertemos o 'IntegerMod_gmp' para int
            
            # Usar power_mod diretamente com 'n' como um número inteiro
            receiver_public.append(power_mod(int(sender_public[choice]), int(self.receiver_values[i]), int(n)))
        
        return receiver_public
    
    def sender_encrypt(self, receiver_public):
        """
        Sender encripta as mensagens
        
        Args:
            receiver_public: Valores enviados pelo receiver
        
        Returns:
            Mensagens encriptadas
        """
        encrypted_messages = []
        
        # Iterar sobre todas as possíveis combinações de duas escolhas
        all_pairs = [(0,1), (0,2), (1,2)]
        
        for i, pair in enumerate(all_pairs):
            # Verificar se este par corresponde às escolhas do receiver
            if pair == tuple(self.choices):
                # Calcular chaves compartilhadas
                shared_key1 = power_mod(int(receiver_public[0]), int(self.sender_keys[pair[0]]), int(self.field.order()))
                shared_key2 = power_mod(int(receiver_public[1]), int(self.sender_keys[pair[1]]), int(self.field.order()))
                
                # Encriptar mensagens usando XOR (simulação simplificada)
                # Em uma implementação real, usaríamos KDF e uma cifra simétrica
                key_hash1 = int(hashlib.sha256(str(shared_key1).encode()).hexdigest(), 16) % 2^32
                key_hash2 = int(hashlib.sha256(str(shared_key2).encode()).hexdigest(), 16) % 2^32
                
                # Codificar as mensagens (na prática, seria uma serialização adequada)
                encoded_msg1 = self.encode_message(self.messages[pair[0]])
                encoded_msg2 = self.encode_message(self.messages[pair[1]])
                
                # Encriptar usando XOR
                encrypted_msg1 = encoded_msg1 ^ key_hash1
                encrypted_msg2 = encoded_msg2 ^ key_hash2
                
                encrypted_messages.append((encrypted_msg1, encrypted_msg2))
            else:
                # Para outras combinações, enviar valores aleatórios
                encrypted_messages.append((secrets.randbits(32), secrets.randbits(32)))
        
        # Embaralhar a ordem se necessário
        # No nosso caso, mantemos a ordem para simplificar
        
        return encrypted_messages
    
    def receiver_decrypt(self, encrypted_messages):
        """
        Receiver decripta as mensagens escolhidas
        
        Args:
            encrypted_messages: Mensagens encriptadas enviadas pelo sender
        
        Returns:
            Mensagens decriptadas
        """
        # Determinar qual índice na lista de pares corresponde às escolhas do receiver
        all_pairs = [(0,1), (0,2), (1,2)]
        pair_idx = all_pairs.index(tuple(self.choices))
        
        # Obter as mensagens encriptadas correspondentes
        encrypted_pair = encrypted_messages[pair_idx]
        
        # Recalcular as chaves compartilhadas
        shared_key1 = power_mod(int(self.g), int(self.sender_keys[self.choices[0]]) * int(self.receiver_values[0]), int(self.field.order()))
        shared_key2 = power_mod(int(self.g), int(self.sender_keys[self.choices[1]]) * int(self.receiver_values[1]), int(self.field.order()))
        
        # Derivar chaves de criptografia
        key_hash1 = int(hashlib.sha256(str(shared_key1).encode()).hexdigest(), 16) % 2^32
        key_hash2 = int(hashlib.sha256(str(shared_key2).encode()).hexdigest(), 16) % 2^32
        
        # Decriptar usando XOR
        decrypted_encoded1 = encrypted_pair[0] ^ key_hash1
        decrypted_encoded2 = encrypted_pair[1] ^ key_hash2
        
        # Decodificar as mensagens
        decrypted_msg1 = self.decode_message(decrypted_encoded1)
        decrypted_msg2 = self.decode_message(decrypted_encoded2)
        
        return [decrypted_msg1, decrypted_msg2]
    
    def encode_message(self, message):
        """
        Codifica uma mensagem (ViewShare) em uma representação numérica
        Simplificação para o exemplo - em uma implementação real seria mais complexo
        """
        # Serializar a mensagem (aqui apenas usando um hash dos componentes)
        # Em uma implementação real, usaríamos uma serialização adequada
        w_hash = sum(int(x) << i for i, x in enumerate(message.w_share))
        z_hash = sum(int(x) << i for i, x in enumerate(message.z_share))
        m_hash = sum(int(x) << i for i, x in enumerate(message.m_share))
        
        return (w_hash * 2^16) + (z_hash * 2^8) + m_hash
    
    def decode_message(self, encoded):
        """
        Decodifica uma representação numérica em uma mensagem (ViewShare)
        Simplificação para o exemplo - em uma implementação real seria mais complexo
        """
        # Extrair componentes (simplificação)
        m_hash = encoded & 0xFF
        z_hash = (encoded >> 8) & 0xFF
        w_hash = (encoded >> 16) & 0xFFFF
        
        # Reconstruir a ViewShare (simplificação - na prática seria mais complexo)
        w_share = [int((w_hash >> i) & 1) for i in range(8)]  # Assumindo 8 bits
        z_share = [int((z_hash >> i) & 1) for i in range(8)]  # Assumindo 8 bits
        m_share = [int((m_hash >> i) & 1) for i in range(8)]  # Assumindo 8 bits
        
        return ViewShare(w_share, z_share, m_share)   

# ======= Parte 5: Utilização do Protocolo =======
class ZKProtocol:
    def __init__(self, circuit, security_param=256, rounds=40):
        self.circuit = circuit
        self.security_param = security_param
        self.rounds = rounds  # Número de rodadas para diminuir o erro de correção
        self.prover = MPCitHProver(circuit, security_param)
        self.ot = ObliviousTransfer(security_param)
    
    def prove(self, witness):
        """
        Executa o protocolo ZK para provar conhecimento do witness
        
        Args:
            witness: Valores de entrada que satisfazem o circuito
        
        Returns:
            Dados da prova
        """
        all_proof_data = []
        
        # Configurar o OT
        p, g = self.ot.setup()
        
        # Executar várias rodadas para diminuir o erro de correção
        for i in range(self.rounds):
            # Gerar um novo sid para esta rodada
            sid = f"session_{i}"
            
            # Prover gera as views dos participantes
            views = self.prover.prove(witness, sid)
            
            # Configurar o OT para esta rodada
            sender_public = self.ot.sender_initial(views)
            
            # Aqui teríamos a interação real com o verificador
            # Para simular, vamos escolher aleatoriamente dois índices
            choices = sorted(random.sample(range(3), 2))
            receiver_public = self.ot.receiver_choose(sender_public, choices)
            
            # Sender encripta as mensagens
            encrypted_messages = self.ot.sender_encrypt(receiver_public)
            
            # Receiver decripta as mensagens
            decrypted_views = self.ot.receiver_decrypt(encrypted_messages)
            
            # Armazenar dados da prova para esta rodada
            proof_data = {
                'sid': sid,
                'sender_public': sender_public,
                'receiver_public': receiver_public,
                'encrypted_messages': encrypted_messages,
                'choices': choices,
                'decrypted_views': decrypted_views
            }
            all_proof_data.append(proof_data)
        
        return all_proof_data
    
    def verify(self, proof_data):
        """
        Verifica se a prova é válida
        
        Args:
            proof_data: Dados da prova gerados pelo prover
        
        Returns:
            True se a prova for válida, False caso contrário
        """
        for round_data in proof_data:
            # Extrair dados da rodada
            decrypted_views = round_data['decrypted_views']
            
            # Verificar consistência entre as views
            if not self.check_consistency(decrypted_views):
                return False
            
            # Verificar se o output do circuito é 1
            if not self.check_output(decrypted_views):
                return False
        
        # Se todas as rodadas forem válidas, aceitar a prova
        return True
    
    def check_consistency(self, views):
        """
        Verifica se as views são consistentes entre si
        
        Args:
            views: Duas views decriptadas
        
        Returns:
            True se as views forem consistentes, False caso contrário
        """
        # Implementação simplificada - verificar se os inputs compartilhados são consistentes
        # Em uma implementação real, verificaríamos todos os wirevalues
        if len(views[0].w_share) != len(views[1].w_share):
            return False
        
        # Para cada par de shares, verificar se são consistentes
        for i in range(len(views[0].w_share)):
            # Implementação simplificada da verificação
            # Em uma implementação real, usaríamos o algoritmo de reconstrução
            pass
            
        return True
    
    def check_output(self, views):
        """
        Verifica se o output do circuito é 1
        
        Args:
            views: Duas views decriptadas
        
        Returns:
            True se o output for 1, False caso contrário
        """
        # Obter o índice do fio de saída
        output_wire = self.circuit.get_output_wire()
        
        # Obter os shares deste fio para as duas views
        output_share1 = views[0].w_share[output_wire] if output_wire < len(views[0].w_share) else 0
        output_share2 = views[1].w_share[output_wire] if output_wire < len(views[1].w_share) else 0
        
        # Em um esquema 2-out-of-3, a soma dos shares deve ser igual ao valor real
        # Simplificação - na prática, usaríamos o algoritmo de reconstrução
        reconstructed_output = (output_share1 + output_share2) % 2
        
        # Verifica se o output é 1
        return reconstructed_output == 1
    
# ======= Parte 7: Main =======
if __name__ == "__main__":
    # Configurações
    lambda_param = 256  # tamanho da seed em bits
    n = 8  # número de entradas do circuito
    
    # Gerar uma seed aleatória
    seed = Integer(secrets.randbits(lambda_param))
    print(f"Seed gerada: {seed}")
    
    # Gerar o circuito
    circuit = generate_polynomial_circuit(seed, n)
    print(circuit)

    # Tenta encontrar um witness aleatório que satisfaça o circuito (output == 1)
    max_attempts = 10000
    for _ in range(max_attempts):
        witness = [randint(0, 1) for _ in range(n)]
        if evaluate_circuit(circuit, witness) == 1:
            break
    else:
        raise ValueError("Não foi possível encontrar um witness válido após várias tentativas.")

    print(f"Witness gerado: {witness}")
    print(f"Avaliação do circuito: {evaluate_circuit(circuit, witness)}")
    
    # Executar o protocolo ZK
    zk = ZKProtocol(circuit, lambda_param, rounds=3)  # Usar apenas 3 rodadas para exemplo
    
    # Gerar prova
    proof_data = zk.prove(witness)
    print(f"Prova gerada com {len(proof_data)} rodadas")
    
    # Verificar prova
    is_valid = zk.verify(proof_data)