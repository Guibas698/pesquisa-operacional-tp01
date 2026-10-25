#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
UFC - Campus Quixadá
QXD0181 - Pesquisa Operacional - Trabalho Prático 01
Prof. Paulo Henrique Macêdo

FORMA 1 (Versão Executável no Windows sem precisar de GCC/CPLEX local):
Abordagem 1 de implementação do algoritmo principal:
  - Na Atividade 1: Validação por Grafo de Adjacência e Busca em Profundidade (DFS)
  - Nas Atividades 2 e 3: Construção Orientada a Objetos (por dicionário de
    variáveis nomeadas x[i,j] e u[i], espelhando linha por linha o modelo
    algébrico do CPLEX / PDF) e resolução via HiGHS.
================================================================================
"""

import os
import sys
import time
import glob
from typing import List, Tuple, Dict, Optional
import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix


class ModeloAlgebricoCPLEXStyle:
    """
    Emula a construção algébrica do Concert Technology (IloModel / IloNumVar)
    em Python puro para que o modelo possa ser lido exatamente como no PDF
    e executado com 1 clique ("Play") no VS Code sem precisar do GCC/CPLEX.
    """

    def __init__(self):
        self.var_indices: Dict[str, int] = {}
        self.lb: List[float] = []
        self.ub: List[float] = []
        self.obj: List[float] = []
        self.is_int: List[int] = []
        self.restricoes: List[Tuple[Dict[int, float], float, float, str]] = []

    def criar_var(
        self, nome: str, lb: float, ub: float, custo_obj: float = 0.0, inteira: bool = False
    ) -> int:
        idx = len(self.lb)
        self.var_indices[nome] = idx
        self.lb.append(float(lb))
        self.ub.append(float(ub))
        self.obj.append(float(custo_obj))
        self.is_int.append(1 if inteira else 0)
        return idx

    def adicionar_igualdade(self, termos: Dict[int, float], valor: float, nome: str = ""):
        self.restricoes.append((termos, float(valor), float(valor), nome))

    def adicionar_menor_igual(self, termos: Dict[int, float], rhs: float, nome: str = ""):
        self.restricoes.append((termos, -np.inf, float(rhs), nome))

    def compilar_matrizes(self):
        num_vars = len(self.lb)
        num_rest = len(self.restricoes)
        A = lil_matrix((num_rest, num_vars), dtype=float)
        lhs = np.zeros(num_rest, dtype=float)
        rhs = np.zeros(num_rest, dtype=float)

        for r, (termos, l_val, r_val, _) in enumerate(self.restricoes):
            for col_idx, coef in termos.items():
                A[r, col_idx] = coef
            lhs[r] = l_val
            rhs[r] = r_val

        return (
            np.array(self.obj, dtype=float),
            A.tocsc(),
            lhs,
            rhs,
            np.array(self.lb, dtype=float),
            np.array(self.ub, dtype=float),
            np.array(self.is_int, dtype=int),
        )


def ler_atsp(caminho: str) -> Tuple[int, List[List[float]]]:
    """Lê a dimensão n e a matriz de custos c[i][j] de um arquivo .atsp."""
    n = 0
    lendo = False
    nums: List[float] = []
    with open(caminho, "r", encoding="utf-8", errors="ignore") as f:
        for linha in f:
            s = linha.strip()
            if not s or s.upper() == "EOF":
                if s.upper() == "EOF":
                    break
                continue
            if not lendo:
                if s.upper().startswith("DIMENSION"):
                    n = int(s.split(":")[-1].strip())
                elif s.upper().startswith("EDGE_WEIGHT_SECTION"):
                    lendo = True
            else:
                for tok in s.split():
                    if tok.upper() == "EOF":
                        break
                    nums.append(float(tok))
    c = [[nums[i * n + j] for j in range(n)] for i in range(n)]
    return n, c


# ==============================================================================
# ATIVIDADE 1 (Forma 1: Busca em Profundidade / Grafo de Adjacência)
# ==============================================================================
def verificar_caminho_valido_dfs(
    arcos: List[Tuple[int, int]],
    n: int,
    s: int,
    t: int,
    c: Optional[List[List[float]]] = None,
) -> Tuple[bool, float]:
    """
    Verifica se a lista de arcos forma um Caminho Hamiltoniano válido de s até t
    construindo uma lista de adjacência e realizando uma travessia única de s a t.
    """
    if len(arcos) != n - 1 or s == t:
        return False, 0.0

    adj: Dict[int, List[int]] = {v: [] for v in range(1, n + 1)}
    grau_in: Dict[int, int] = {v: 0 for v in range(1, n + 1)}
    custo = 0.0

    for u, v in arcos:
        if u < 1 or u > n or v < 1 or v > n or u == v:
            return False, 0.0
        adj[u].append(v)
        grau_in[v] += 1
        if c is not None:
            custo += c[u - 1][v - 1]

    # Graus obrigatórios (Equações 3 a 8)
    if len(adj[s]) != 1 or grau_in[s] != 0:
        return False, 0.0
    if grau_in[t] != 1 or len(adj[t]) != 0:
        return False, 0.0
    for k in range(1, n + 1):
        if k not in (s, t) and (len(adj[k]) != 1 or grau_in[k] != 1):
            return False, 0.0

    # Travessia DFS a partir da origem s
    visitados = set()
    atual = s
    while True:
        if atual in visitados:
            return False, 0.0
        visitados.add(atual)
        if atual == t:
            break
        atual = adj[atual][0]

    return (len(visitados) == n and atual == t), custo


# ==============================================================================
# ATIVIDADES 2 E 3 (Forma 1: Construção Algébrica estilo CPLEX)
# ==============================================================================
def montar_e_resolver_mtz(
    n: int, c: List[List[float]], relaxacao_linear: bool, limite_tempo: float = 30.0
):
    s, t = 0, n - 1
    modelo = ModeloAlgebricoCPLEXStyle()

    # 1. Criação das variáveis x_{i,j} (Eq. 1)
    x = {}
    for i in range(n):
        for j in range(n):
            if i == j:
                x[i, j] = modelo.criar_var(f"x_{i+1}_{j+1}", 0.0, 0.0, 0.0, inteira=False)
            else:
                x[i, j] = modelo.criar_var(
                    f"x_{i+1}_{j+1}",
                    0.0,
                    1.0,
                    custo_obj=c[i][j],
                    inteira=(not relaxacao_linear),
                )

    # 2. Criação das variáveis de ordem de visitação u_i (Eq. 9 e Eq. 11)
    u = {}
    for i in range(n):
        if i == s:
            u[i] = modelo.criar_var(f"u_{i+1}", 1.0, 1.0, 0.0, inteira=False)
        else:
            u[i] = modelo.criar_var(f"u_{i+1}", 2.0, float(n), 0.0, inteira=False)

    # Eq. (3): O caminho deve sair de s
    modelo.adicionar_igualdade({x[s, j]: 1.0 for j in range(n) if j != s}, 1.0, "saida_s")

    # Eq. (4): O caminho deve chegar em t
    modelo.adicionar_igualdade({x[i, t]: 1.0 for i in range(n) if i != t}, 1.0, "chegada_t")

    # Eq. (5): O caminho não pode retornar a s
    modelo.adicionar_igualdade({x[i, s]: 1.0 for i in range(n) if i != s}, 0.0, "nao_retorna_s")

    # Eq. (6): O caminho não pode sair de t
    modelo.adicionar_igualdade({x[t, j]: 1.0 for j in range(n) if j != t}, 0.0, "nao_sai_t")

    # Eq. (7) e Eq. (8): Passagem única pelos vértices intermediários k in V \ {s, t}
    for k in range(n):
        if k in (s, t):
            continue
        modelo.adicionar_igualdade({x[i, k]: 1.0 for i in range(n) if i != k}, 1.0, f"in_{k+1}")
        modelo.adicionar_igualdade({x[k, j]: 1.0 for j in range(n) if j != k}, 1.0, f"out_{k+1}")

    # Eq. (10): Restrições MTZ: u_i - u_j + n * x_{i,j} <= n - 1
    for i in range(n):
        for j in range(n):
            if j == s or i == j:
                continue
            modelo.adicionar_menor_igual(
                {u[i]: 1.0, u[j]: -1.0, x[i, j]: float(n)},
                float(n - 1),
                f"mtz_{i+1}_{j+1}",
            )

    obj, A_csc, lhs, rhs, lb, ub, integrality = modelo.compilar_matrizes()

    inicio = time.perf_counter()
    if relaxacao_linear:
        mask_eq = np.isclose(lhs, rhs)
        res = linprog(
            c=obj,
            A_ub=A_csc[~mask_eq, :],
            b_ub=rhs[~mask_eq],
            A_eq=A_csc[mask_eq, :],
            b_eq=rhs[mask_eq],
            bounds=list(zip(lb, ub)),
            method="highs",
        )
        tempo = time.perf_counter() - inicio
        return res.fun, tempo, []
    else:
        res = milp(
            c=obj,
            integrality=integrality,
            bounds=Bounds(lb, ub),
            constraints=LinearConstraint(A_csc, lhs, rhs),
            options={"time_limit": limite_tempo},
        )
        tempo = time.perf_counter() - inicio
        arcos_sol = []
        if res.x is not None:
            prox = {}
            for i in range(n):
                for j in range(n):
                    if i != j and res.x[x[i, j]] > 0.5:
                        prox[i + 1] = j + 1
            atual = 1
            while atual != n and len(arcos_sol) < n:
                seg = prox.get(atual)
                if seg is None:
                    break
                arcos_sol.append((atual, seg))
                atual = seg
        return res.fun, tempo, arcos_sol


def main():
    pasta_atual = os.path.dirname(os.path.abspath(__file__))
    pasta_inst = os.path.normpath(
        os.path.join(
            pasta_atual,
            "..",
            "..",
            "Implementacao",
            "instancias",
            "Completos e Assimetricos nos Arcos",
        )
    )

    arquivos = sorted(glob.glob(os.path.join(pasta_inst, "*.atsp")))
    if len(sys.argv) > 1 and not sys.argv[1].startswith("--"):
        arquivos = [sys.argv[1]]

    print("=" * 92)
    print("EXECUÇÃO FORMA 1 (Modelagem Algébrica Estilo CPLEX + Validação por DFS)")
    print("=" * 92)
    print(f"{'Instância':<12} | {'n':>3} | {'Relax. Linear':>14} | {'Sol. Inteira':>14} | {'Tempo PLI':>9} | {'Validação'}")
    print("-" * 92)

    for arq in arquivos:
        nome = os.path.basename(arq)
        n, c = ler_atsp(arq)
        val_lp, _, _ = montar_e_resolver_mtz(n, c, relaxacao_linear=True)
        val_pli, t_pli, arcos = montar_e_resolver_mtz(n, c, relaxacao_linear=False, limite_tempo=30.0)
        ok, _ = verificar_caminho_valido_dfs(arcos, n, 1, n, c)
        print(
            f"{nome:<12} | {n:>3} | {val_lp:>14.4f} | {val_pli:>14.2f} | {t_pli:>8.4f}s | {'VÁLIDO' if ok else 'INVÁLIDO'}"
        )
    print("=" * 92)


if __name__ == "__main__":
    main()
