#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
UFC - Campus Quixadá
QXD0181 - Pesquisa Operacional - Trabalho Prático 01
Prof. Paulo Henrique Macêdo

Programação Linear Inteira para o Problema do Caminho Hamiltoniano Mínimo (PCHM)
com Vértices de Origem (s = 1) e Destino (t = n) - Formulação MTZ

Implementação Forma 2: Python (com solver HiGHS via SciPy / PuLP)
Contempla integralmente as Atividades 1, 2 e 3 descritas no PDF e executa
automaticamente todas as 13 instâncias disponibilizadas.
================================================================================
"""

import os
import sys
import time
import glob
import csv
from typing import List, Tuple, Optional, Dict, Any
import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix


# ==============================================================================
# LEITOR DE INSTÂNCIAS TSPLIB (.atsp - FULL_MATRIX)
# ==============================================================================
def ler_instancia_atsp(caminho_arquivo: str) -> Tuple[str, int, np.ndarray]:
    """
    Lê um arquivo .atsp no formato TSPLIB (EDGE_WEIGHT_FORMAT: FULL_MATRIX)
    e retorna (nome, n, matriz_custos_n_por_n).
    """
    nome = os.path.basename(caminho_arquivo)
    n = 0
    lendo_pesos = False
    valores: List[float] = []

    with open(caminho_arquivo, "r", encoding="utf-8", errors="ignore") as f:
        for linha in f:
            linha_limpa = linha.strip()
            if not linha_limpa:
                continue
            if linha_limpa.upper() == "EOF":
                break

            if not lendo_pesos:
                if linha_limpa.upper().startswith("NAME"):
                    partes = linha_limpa.split(":", 1)
                    if len(partes) > 1:
                        nome = partes[1].strip()
                elif linha_limpa.upper().startswith("DIMENSION"):
                    partes = linha_limpa.split(":", 1)
                    if len(partes) > 1:
                        n = int(partes[1].strip())
                elif linha_limpa.upper().startswith("EDGE_WEIGHT_SECTION"):
                    lendo_pesos = True
            else:
                for token in linha_limpa.split():
                    if token.upper() == "EOF":
                        break
                    valores.append(float(token))

    if n <= 0 or len(valores) < n * n:
        raise ValueError(
            f"Erro ao ler instância {caminho_arquivo}: n={n}, lidos={len(valores)} pesos."
        )

    matriz_custos = np.array(valores[: n * n], dtype=float).reshape((n, n))
    return nome, n, matriz_custos


# ==============================================================================
# ATIVIDADE 1:
# "Implemente uma função que receba uma sequência de arcos do grafo e verifique
# se ela corresponde a um caminho válido, sendo uma solução viável para o
# problema do Caminho Hamiltoniano Mínimo (PCHM) com Vértices de Origem e Destino."
# ==============================================================================
def verificar_caminho_valido(
    arcos: List[Tuple[int, int]],
    n: int,
    s: int = 1,
    t: Optional[int] = None,
    matriz_custos: Optional[np.ndarray] = None,
    exibir_detalhes: bool = False,
) -> Tuple[bool, str, float]:
    """
    Verifica se uma sequência de arcos [(u_1, v_1), (u_2, v_2), ...] (com vértices
    numerados de 1 até n) corresponde a um Caminho Hamiltoniano válido que inicia
    em s (padrão 1), termina em t (padrão n) e visita todos os vértices de
    V = {1, ..., n} exatamente uma vez.

    Retorna:
        (eh_valido: bool, mensagem: str, custo_total: float)
    """
    if t is None:
        t = n

    # 1. Validação dos parâmetros básicos do problema
    if n < 2:
        msg = f"Número de vértices inválido (n={n})."
        if exibir_detalhes:
            print(msg)
        return False, msg, 0.0

    if not (1 <= s <= n) or not (1 <= t <= n) or s == t:
        msg = f"Vértices de origem (s={s}) e destino (t={t}) inválidos ou iguais."
        if exibir_detalhes:
            print(msg)
        return False, msg, 0.0

    # 2. Quantidade de arcos deve ser exatamente n - 1
    if len(arcos) != n - 1:
        msg = f"Quantidade de arcos ({len(arcos)}) diferente de n - 1 ({n - 1})."
        if exibir_detalhes:
            print(msg)
        return False, msg, 0.0

    grau_entrada = {v: 0 for v in range(1, n + 1)}
    grau_saida = {v: 0 for v in range(1, n + 1)}
    proximo_vertice: Dict[int, int] = {}
    custo_total = 0.0

    # 3. Verifica cada arco individualmente e verifica continuidade da sequência
    for idx, arco in enumerate(arcos):
        if len(arco) != 2:
            msg = f"Elemento na posição {idx} não é um par (u, v): {arco}."
            if exibir_detalhes:
                print(msg)
            return False, msg, 0.0

        u, v = arco
        if not (1 <= u <= n) or not (1 <= v <= n):
            msg = f"Arco ({u}, {v}) contém vértice fora do conjunto V={{1, ..., {n}}}."
            if exibir_detalhes:
                print(msg)
            return False, msg, 0.0

        if u == v:
            msg = f"Auto-laço detectado no vértice {u}: arco ({u}, {v})."
            if exibir_detalhes:
                print(msg)
            return False, msg, 0.0

        grau_saida[u] += 1
        grau_entrada[v] += 1
        proximo_vertice[u] = v

        if matriz_custos is not None:
            custo_total += float(matriz_custos[u - 1, v - 1])

    # 4. Verifica se a sequência passada já está encadeada ( (v_0,v_1), (v_1,v_2), ... )
    if arcos[0][0] != s:
        msg = f"A sequência de arcos não inicia no vértice de origem s={s} (iniciou em {arcos[0][0]})."
        if exibir_detalhes:
            print(msg)
        return False, msg, 0.0

    if arcos[-1][1] != t:
        msg = f"A sequência de arcos não termina no vértice de destino t={t} (terminou em {arcos[-1][1]})."
        if exibir_detalhes:
            print(msg)
        return False, msg, 0.0

    for idx in range(len(arcos) - 1):
        if arcos[idx][1] != arcos[idx + 1][0]:
            msg = (
                f"Descontinuidade na sequência entre o arco {idx} {arcos[idx]} "
                f"e o arco {idx + 1} {arcos[idx + 1]}."
            )
            if exibir_detalhes:
                print(msg)
            return False, msg, 0.0

    # 5. Verifica as restrições de grau (Equações 3 a 8 do PDF)
    if grau_saida[s] != 1 or grau_entrada[s] != 0:
        msg = f"Violação de grau na origem s={s}: saída={grau_saida[s]} (esperado 1), entrada={grau_entrada[s]} (esperado 0)."
        if exibir_detalhes:
            print(msg)
        return False, msg, 0.0

    if grau_entrada[t] != 1 or grau_saida[t] != 0:
        msg = f"Violação de grau no destino t={t}: entrada={grau_entrada[t]} (esperado 1), saída={grau_saida[t]} (esperado 0)."
        if exibir_detalhes:
            print(msg)
        return False, msg, 0.0

    for k in range(1, n + 1):
        if k in (s, t):
            continue
        if grau_entrada[k] != 1 or grau_saida[k] != 1:
            msg = (
                f"Violação de conservação de fluxo no vértice intermediário k={k}: "
                f"entrada={grau_entrada[k]} (esperado 1), saída={grau_saida[k]} (esperado 1)."
            )
            if exibir_detalhes:
                print(msg)
            return False, msg, 0.0

    # 6. Travessia do caminho de s até t para garantir ausência de subrotas / ciclos
    visitados = set()
    atual = s
    while atual is not None:
        if atual in visitados:
            msg = f"Ciclo fechado detectado ao revisitar o vértice {atual}."
            if exibir_detalhes:
                print(msg)
            return False, msg, 0.0
        visitados.add(atual)
        if atual == t:
            break
        atual = proximo_vertice.get(atual, None)

    if len(visitados) != n or atual != t:
        msg = f"O caminho não visitou todos os {n} vértices (visitou {len(visitados)})."
        if exibir_detalhes:
            print(msg)
        return False, msg, 0.0

    msg = f"Caminho Hamiltoniano válido de s={s} até t={t} visitando todos os {n} vértices!"
    if exibir_detalhes:
        print(msg)
    return True, msg, custo_total


# ==============================================================================
# CONSTRUÇÃO MATRICIAL DO MODELO MTZ (ATIVIDADES 2 E 3)
# ==============================================================================
def construir_modelo_mtz(n: int, c: np.ndarray, s: int = 1, t: Optional[int] = None):
    """
    Constrói o modelo matemático de Programação Linear Inteira (MTZ) conforme
    as Equações (1) a (11) do PDF.

    Indexação das variáveis no vetor de decisão y de dimensão (n*n + n):
      - Variáveis x_{i,j} para i, j in {0, ..., n-1}: índice = i * n + j
      - Variáveis u_i para i in {0, ..., n-1}: índice = n*n + i
    """
    if t is None:
        t = n
    s_idx = s - 1
    t_idx = t - 1

    num_x = n * n
    num_u = n
    num_vars = num_x + num_u

    # Função Objetivo (Eq. 2): Minimizar sum_{i in V} sum_{j in V} c_{ij} x_{ij}
    obj = np.zeros(num_vars, dtype=float)
    for i in range(n):
        for j in range(n):
            if i != j:
                obj[i * n + j] = c[i, j]

    # Limites das variáveis (Bounds):
    # - x_{ij} in [0, 1] para i != j, e x_{ii} = 0 (sem auto-laços)
    # - Eq. (9): u_s = 1
    # - Eq. (11): 2 <= u_i <= n, para todo i in V \ {s}
    lb = np.zeros(num_vars, dtype=float)
    ub = np.ones(num_vars, dtype=float)

    for i in range(n):
        # Auto-laço x_{i,i} fixado em 0
        lb[i * n + i] = 0.0
        ub[i * n + i] = 0.0

        # Limites de u_i
        u_idx = num_x + i
        if i == s_idx:
            lb[u_idx] = 1.0
            ub[u_idx] = 1.0
        else:
            lb[u_idx] = 2.0
            ub[u_idx] = float(n)

    # Lista de restrições: (indices, coeficientes, limite_inferior, limite_superior)
    linhas_idx = []
    linhas_val = []
    lhs_list = []
    rhs_list = []

    def adicionar_restricao(indices, valores, lb_val, ub_val):
        linhas_idx.append(indices)
        linhas_val.append(valores)
        lhs_list.append(lb_val)
        rhs_list.append(ub_val)

    # Eq. (3): sum_{j in V} x_{s,j} = 1 (O caminho deve sair de s)
    idxs = [s_idx * n + j for j in range(n) if j != s_idx]
    adicionar_restricao(idxs, [1.0] * len(idxs), 1.0, 1.0)

    # Eq. (4): sum_{i in V} x_{i,t} = 1 (O caminho deve chegar em t)
    idxs = [i * n + t_idx for i in range(n) if i != t_idx]
    adicionar_restricao(idxs, [1.0] * len(idxs), 1.0, 1.0)

    # Eq. (5): sum_{i in V} x_{i,s} = 0 (O caminho não pode retornar a s)
    idxs = [i * n + s_idx for i in range(n) if i != s_idx]
    adicionar_restricao(idxs, [1.0] * len(idxs), 0.0, 0.0)

    # Eq. (6): sum_{j in V} x_{t,j} = 0 (O caminho não pode sair de t)
    idxs = [t_idx * n + j for j in range(n) if j != t_idx]
    adicionar_restricao(idxs, [1.0] * len(idxs), 0.0, 0.0)

    # Eq. (7) e Eq. (8): Conservação de fluxo nos vértices intermediários k in V \ {s, t}
    for k in range(n):
        if k == s_idx or k == t_idx:
            continue
        # Eq. (7): sum_{i in V} x_{i,k} = 1
        idxs_in = [i * n + k for i in range(n) if i != k]
        adicionar_restricao(idxs_in, [1.0] * len(idxs_in), 1.0, 1.0)

        # Eq. (8): sum_{j in V} x_{k,j} = 1
        idxs_out = [k * n + j for j in range(n) if j != k]
        adicionar_restricao(idxs_out, [1.0] * len(idxs_out), 1.0, 1.0)

    # Eq. (10): Restrições MTZ: u_i - u_j + n * x_{i,j} <= n - 1
    # Para todo i in V, para todo j in V \ {s}, com i != j
    for i in range(n):
        for j in range(n):
            if j == s_idx or i == j:
                continue
            u_i = num_x + i
            u_j = num_x + j
            x_ij = i * n + j
            adicionar_restricao([u_i, u_j, x_ij], [1.0, -1.0, float(n)], -np.inf, float(n - 1))

    num_restricoes = len(lhs_list)
    A = lil_matrix((num_restricoes, num_vars), dtype=float)
    for r, (idxs, vals) in enumerate(zip(linhas_idx, linhas_val)):
        A[r, idxs] = vals

    A_csc = A.tocsc()
    lhs = np.array(lhs_list, dtype=float)
    rhs = np.array(rhs_list, dtype=float)

    return obj, A_csc, lhs, rhs, lb, ub


# ==============================================================================
# ATIVIDADE 2:
# "Implemente o modelo de Programação Linear Inteira com as restrições de
# Miller-Tucker-Zemlin para o PCHM com Vértices de Origem e Destino e resolva
# sua relaxação linear para cada instância disponibilizada."
# ==============================================================================
def resolver_relaxacao_linear(n: int, c: np.ndarray, s: int = 1, t: Optional[int] = None) -> Dict[str, Any]:
    """
    Resolve a relaxação linear (LP) do modelo MTZ para o PCHM.
    """
    if t is None:
        t = n
    obj, A_csc, lhs, rhs, lb, ub = construir_modelo_mtz(n, c, s, t)

    # Para linprog(method='highs'), transformamos lhs <= A x <= rhs em desigualdades A_ub x <= b_ub e igualdades A_eq x = b_eq
    mask_eq = np.isclose(lhs, rhs)
    mask_ub = ~mask_eq

    A_eq = A_csc[mask_eq, :]
    b_eq = rhs[mask_eq]
    A_ub = A_csc[mask_ub, :]
    b_ub = rhs[mask_ub]
    bounds = list(zip(lb, ub))

    inicio = time.perf_counter()
    res = linprog(
        c=obj,
        A_ub=A_ub,
        b_ub=b_ub,
        A_eq=A_eq,
        b_eq=b_eq,
        bounds=bounds,
        method="highs",
    )
    tempo = time.perf_counter() - inicio

    return {
        "sucesso": bool(res.success),
        "status": res.message,
        "valor_objetivo": float(res.fun) if res.success else None,
        "tempo_segundos": tempo,
        "solucao_x": res.x[: n * n].reshape((n, n)) if res.success else None,
        "solucao_u": res.x[n * n :] if res.success else None,
    }


# ==============================================================================
# ATIVIDADE 3:
# "Adicione a restrição de integralidade às variáveis e execute o modelo para
# que retorne a melhor solução inteira encontrada em até 3 minutos de execução
# para cada instância."
# ==============================================================================
def resolver_modelo_inteiro(
    n: int,
    c: np.ndarray,
    s: int = 1,
    t: Optional[int] = None,
    limite_tempo_segundos: float = 180.0,
) -> Dict[str, Any]:
    """
    Resolve o modelo de Programação Linear Inteira (PLI) com restrições MTZ
    e limite de tempo de 3 minutos (180 segundos).
    Extrai o caminho encontrado e valida com a função da Atividade 1.
    """
    if t is None:
        t = n
    obj, A_csc, lhs, rhs, lb, ub = construir_modelo_mtz(n, c, s, t)

    # Variáveis x_{ij} binárias (integrality = 1).
    # Nota matemática importante sobre MTZ (Eq. 9-11):
    # Como u_s = 1, 2 <= u_i <= n e u_j >= u_i + 1 para todo arco (i,j) com x_{ij} = 1,
    # ao longo dos n - 1 arcos do caminho hamiltoniano s = v_1 -> v_2 -> ... -> v_n = t,
    # tem-se 1 = u_{v_1} <= u_{v_2} - 1 <= ... <= u_{v_n} - (n - 1) <= n - (n - 1) = 1.
    # Logo, quando x_{ij} in {0, 1}, u_{v_k} = k (inteiro exato!) é forçado automaticamente.
    # Focar o branch-and-bound nas variáveis binárias x_{ij} evita ramificações redundantes em u_i.
    integrality = np.ones_like(obj)
    integrality[n * n :] = 0
    constraints = LinearConstraint(A_csc, lhs, rhs)
    bounds = Bounds(lb, ub)

    inicio = time.perf_counter()
    res = milp(
        c=obj,
        integrality=integrality,
        bounds=bounds,
        constraints=constraints,
        options={"time_limit": limite_tempo_segundos},
    )
    tempo = time.perf_counter() - inicio

    tem_solucao = res.x is not None and len(res.x) > 0
    if not tem_solucao:
        return {
            "sucesso": False,
            "status": res.message,
            "valor_objetivo": None,
            "gap_percentual": None,
            "tempo_segundos": tempo,
            "arcos": [],
            "caminho_valido": False,
            "custo_verificado": None,
        }

    x_sol = res.x[: n * n].reshape((n, n))
    u_sol = res.x[n * n :]

    # Reconstrói a sequência ordenada de arcos a partir de s
    proximo = {}
    for i in range(n):
        for j in range(n):
            if i != j and x_sol[i, j] > 0.5:
                proximo[i + 1] = j + 1

    arcos_caminho: List[Tuple[int, int]] = []
    atual = s
    while atual != t and len(arcos_caminho) < n:
        seguinte = proximo.get(atual, None)
        if seguinte is None:
            break
        arcos_caminho.append((atual, seguinte))
        atual = seguinte

    eh_valido, msg_validacao, custo_verificado = verificar_caminho_valido(
        arcos_caminho, n, s=s, t=t, matriz_custos=c, exibir_detalhes=False
    )

    gap_relativo = getattr(res, "mip_gap", 0.0)
    if gap_relativo is None or np.isnan(gap_relativo):
        gap_relativo = 0.0

    status_str = "Ótimo" if res.status == 0 else f"Viável (Limite de Tempo / {res.message})"

    return {
        "sucesso": True,
        "status": status_str,
        "valor_objetivo": float(res.fun),
        "gap_percentual": float(gap_relativo) * 100.0,
        "tempo_segundos": tempo,
        "arcos": arcos_caminho,
        "caminho_valido": eh_valido,
        "mensagem_validacao": msg_validacao,
        "custo_verificado": custo_verificado,
        "u_sol": u_sol,
    }


# ==============================================================================
# TESTES UNITÁRIOS DA ATIVIDADE 1 + EXECUÇÃO EM TODAS AS INSTÂNCIAS
# ==============================================================================
def testar_validador_atividade_1():
    print("=" * 80)
    print("DEMONSTRAÇÃO ATIVIDADE 1: Validador de Caminho Hamiltoniano (s=1, t=n)")
    print("=" * 80)
    n_teste = 5
    # Caminho válido: 1 -> 3 -> 2 -> 4 -> 5
    arcos_validos = [(1, 3), (3, 2), (2, 4), (4, 5)]
    ok, msg, _ = verificar_caminho_valido(arcos_validos, n=n_teste, s=1, t=5)
    print(f"Teste 1 (Caminho válido 1->3->2->4->5): {ok} | {msg}")

    # Caminho inválido com subrota (1 -> 5 direto e ciclo 2 -> 3 -> 4 -> 2)
    arcos_subrota = [(1, 5), (2, 3), (3, 4), (4, 2)]
    ok2, msg2, _ = verificar_caminho_valido(arcos_subrota, n=n_teste, s=1, t=5)
    print(f"Teste 2 (Com subrota disjunta)        : {ok2} | {msg2}")

    # Caminho inválido que não visita todos os vértices (1 -> 2 -> 5)
    arcos_incompleto = [(1, 2), (2, 5)]
    ok3, msg3, _ = verificar_caminho_valido(arcos_incompleto, n=n_teste, s=1, t=5)
    print(f"Teste 3 (Caminho incompleto 1->2->5)  : {ok3} | {msg3}")
    print()


def main():
    testar_validador_atividade_1()

    # Localiza pasta de instâncias
    pasta_base = os.path.dirname(os.path.abspath(__file__))
    pasta_instancias = os.path.normpath(
        os.path.join(
            pasta_base,
            "..",
            "..",
            "Implementacao",
            "instancias",
            "Completos e Assimetricos nos Arcos",
        )
    )

    limite_tempo = 180.0
    args_arquivos = []
    for arg in sys.argv[1:]:
        if arg == "--rapido":
            limite_tempo = 15.0
        elif arg.startswith("--limite="):
            limite_tempo = float(arg.split("=", 1)[1])
        else:
            args_arquivos.append(arg)

    if args_arquivos:
        arquivos = args_arquivos
    else:
        arquivos = sorted(glob.glob(os.path.join(pasta_instancias, "*.atsp")))

    if not arquivos:
        print(f"Nenhuma instância encontrada em: {pasta_instancias}")
        return

    resultados_gerais = []

    print("=" * 105)
    print(
        f"{'Instância':<12} | {'n':>3} | {'Relax. Linear':>14} | {'Tempo LP':>9} | "
        f"{'Sol. Inteira':>14} | {'Gap (%)':>8} | {'Tempo PLI':>9} | {'Validação (Atv.1)'}"
    )
    print("-" * 105)

    for caminho in arquivos:
        nome_arq = os.path.basename(caminho)
        nome, n, c = ler_instancia_atsp(caminho)
        s, t = 1, n

        # Atividade 2: Relaxação Linear
        res_lp = resolver_relaxacao_linear(n, c, s=s, t=t)

        # Atividade 3: Modelo Inteiro (até 3 minutos)
        res_pli = resolver_modelo_inteiro(n, c, s=s, t=t, limite_tempo_segundos=limite_tempo)

        lp_val = res_lp["valor_objetivo"]
        pli_val = res_pli["valor_objetivo"]
        gap_pli = res_pli["gap_percentual"]
        valido_str = "VÁLIDO" if res_pli["caminho_valido"] else "INVÁLIDO"

        lp_str = f"{lp_val:>14.4f}" if lp_val is not None else f"{'N/A':>14}"
        pli_str = f"{pli_val:>14.2f}" if pli_val is not None else f"{'Sem sol.':>14}"
        gap_str = f"{gap_pli:>7.2f}%" if gap_pli is not None else f"{'N/A':>8}"

        print(
            f"{nome_arq:<12} | {n:>3} | {lp_str} | {res_lp['tempo_segundos']:>8.4f}s | "
            f"{pli_str} | {gap_str} | {res_pli['tempo_segundos']:>8.4f}s | {valido_str}"
        )
        sys.stdout.flush()

        caminho_str = " -> ".join([str(s)] + [str(v) for (_, v) in res_pli["arcos"]])

        resultados_gerais.append(
            {
                "instancia": nome_arq,
                "n": n,
                "s": s,
                "t": t,
                "relaxacao_linear": round(lp_val, 4) if lp_val is not None else "",
                "tempo_lp_s": round(res_lp["tempo_segundos"], 4),
                "solucao_inteira": round(pli_val, 2) if pli_val is not None else "",
                "gap_mip_pct": round(gap_pli, 2) if gap_pli is not None else "",
                "tempo_pli_s": round(res_pli["tempo_segundos"], 4),
                "status_pli": res_pli["status"],
                "caminho_valido_atv1": res_pli["caminho_valido"],
                "custo_verificado_atv1": round(res_pli["custo_verificado"], 2)
                if res_pli["custo_verificado"] is not None
                else "",
                "caminho": caminho_str,
            }
        )

    print("=" * 105)

    # Salva CSV com todos os resultados
    caminho_csv = os.path.join(pasta_base, "resultados_instancias.csv")
    with open(caminho_csv, "w", newline="", encoding="utf-8") as f_csv:
        writer = csv.DictWriter(f_csv, fieldnames=list(resultados_gerais[0].keys()))
        writer.writeheader()
        writer.writerows(resultados_gerais)

    print(f"\nResultados salvos com sucesso em: {caminho_csv}")


if __name__ == "__main__":
    main()
