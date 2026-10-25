import os
import sys
import time
import glob
import csv
import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix


def ler_instancia_atsp(caminho_arquivo):
    nome = os.path.basename(caminho_arquivo)
    n = 0
    lendo_pesos = False
    valores = []

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

    matriz_custos = np.array(valores[: n * n], dtype=float).reshape((n, n))
    return nome, n, matriz_custos


# Questao 1: verifica se a sequencia de arcos corresponde a um caminho hamiltoniano valido de s a t
def verificar_caminho_valido(arcos, n, s=1, t=None, matriz_custos=None):
    if t is None:
        t = n

    if n < 2 or not (1 <= s <= n) or not (1 <= t <= n) or s == t:
        return False, "Parametros de origem/destino invalidos.", 0.0

    if len(arcos) != n - 1:
        return False, f"Quantidade de arcos ({len(arcos)}) diferente de n - 1 ({n - 1}).", 0.0

    grau_entrada = {v: 0 for v in range(1, n + 1)}
    grau_saida = {v: 0 for v in range(1, n + 1)}
    proximo_vertice = {}
    custo_total = 0.0

    for u, v in arcos:
        if not (1 <= u <= n) or not (1 <= v <= n) or u == v:
            return False, f"Arco invalido: ({u}, {v}).", 0.0

        grau_saida[u] += 1
        grau_entrada[v] += 1
        proximo_vertice[u] = v

        if matriz_custos is not None:
            custo_total += float(matriz_custos[u - 1, v - 1])

    if arcos[0][0] != s or arcos[-1][1] != t:
        return False, "O caminho nao inicia em s ou nao termina em t.", 0.0

    for idx in range(len(arcos) - 1):
        if arcos[idx][1] != arcos[idx + 1][0]:
            return False, "Sequencia de arcos desconexa.", 0.0

    # Verifica graus (Eq. 3 a 8)
    if grau_saida[s] != 1 or grau_entrada[s] != 0:
        return False, "Grau invalido na origem s.", 0.0
    if grau_entrada[t] != 1 or grau_saida[t] != 0:
        return False, "Grau invalido no destino t.", 0.0

    for k in range(1, n + 1):
        if k in (s, t):
            continue
        if grau_entrada[k] != 1 or grau_saida[k] != 1:
            return False, f"Grau invalido no vertice {k}.", 0.0

    # Percorre de s ate t para checar ciclos/subrotas
    visitados = set()
    atual = s
    while atual is not None:
        if atual in visitados:
            return False, "Ciclo detectado.", 0.0
        visitados.add(atual)
        if atual == t:
            break
        atual = proximo_vertice.get(atual, None)

    if len(visitados) != n or atual != t:
        return False, "Nem todos os vertices foram visitados.", 0.0

    return True, "Caminho valido", custo_total


def construir_modelo_mtz(n, c, s=1, t=None):
    if t is None:
        t = n
    s_idx = s - 1
    t_idx = t - 1

    num_x = n * n
    num_u = n
    num_vars = num_x + num_u

    # Funcao objetivo (Eq. 2)
    obj = np.zeros(num_vars, dtype=float)
    for i in range(n):
        for j in range(n):
            if i != j:
                obj[i * n + j] = c[i, j]

    lb = np.zeros(num_vars, dtype=float)
    ub = np.ones(num_vars, dtype=float)

    for i in range(n):
        lb[i * n + i] = 0.0
        ub[i * n + i] = 0.0

        u_idx = num_x + i
        if i == s_idx:
            lb[u_idx] = 1.0
            ub[u_idx] = 1.0
        else:
            lb[u_idx] = 2.0
            ub[u_idx] = float(n)

    linhas_idx = []
    linhas_val = []
    lhs_list = []
    rhs_list = []

    def adicionar_restricao(indices, valores, lb_val, ub_val):
        linhas_idx.append(indices)
        linhas_val.append(valores)
        lhs_list.append(lb_val)
        rhs_list.append(ub_val)

    # Eq. (3) a (6): graus em s e t
    idxs = [s_idx * n + j for j in range(n) if j != s_idx]
    adicionar_restricao(idxs, [1.0] * len(idxs), 1.0, 1.0)

    idxs = [i * n + t_idx for i in range(n) if i != t_idx]
    adicionar_restricao(idxs, [1.0] * len(idxs), 1.0, 1.0)

    idxs = [i * n + s_idx for i in range(n) if i != s_idx]
    adicionar_restricao(idxs, [1.0] * len(idxs), 0.0, 0.0)

    idxs = [t_idx * n + j for j in range(n) if j != t_idx]
    adicionar_restricao(idxs, [1.0] * len(idxs), 0.0, 0.0)

    # Eq. (7) e (8): graus nos vertices intermediarios
    for k in range(n):
        if k == s_idx or k == t_idx:
            continue
        idxs_in = [i * n + k for i in range(n) if i != k]
        adicionar_restricao(idxs_in, [1.0] * len(idxs_in), 1.0, 1.0)

        idxs_out = [k * n + j for j in range(n) if j != k]
        adicionar_restricao(idxs_out, [1.0] * len(idxs_out), 1.0, 1.0)

    # Eq. (10): restricoes MTZ
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

    return obj, A.tocsc(), np.array(lhs_list, dtype=float), np.array(rhs_list, dtype=float), lb, ub


# Questao 2: resolve a relaxacao linear
def resolver_relaxacao_linear(n, c, s=1, t=None):
    if t is None:
        t = n
    obj, A_csc, lhs, rhs, lb, ub = construir_modelo_mtz(n, c, s, t)

    mask_eq = np.isclose(lhs, rhs)
    mask_ub = ~mask_eq

    inicio = time.perf_counter()
    res = linprog(
        c=obj,
        A_ub=A_csc[mask_ub, :],
        b_ub=rhs[mask_ub],
        A_eq=A_csc[mask_eq, :],
        b_eq=rhs[mask_eq],
        bounds=list(zip(lb, ub)),
        method="highs",
    )
    tempo = time.perf_counter() - inicio

    return {
        "sucesso": bool(res.success),
        "valor_objetivo": float(res.fun) if res.success else None,
        "tempo_segundos": tempo,
    }


# Questao 3: resolve o modelo inteiro com limite de 3 minutos (180s)
def resolver_modelo_inteiro(n, c, s=1, t=None, limite_tempo_segundos=180.0):
    if t is None:
        t = n
    obj, A_csc, lhs, rhs, lb, ub = construir_modelo_mtz(n, c, s, t)

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

    if res.x is None or len(res.x) == 0:
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

    proximo = {}
    for i in range(n):
        for j in range(n):
            if i != j and x_sol[i, j] > 0.5:
                proximo[i + 1] = j + 1

    arcos_caminho = []
    atual = s
    while atual != t and len(arcos_caminho) < n:
        seguinte = proximo.get(atual, None)
        if seguinte is None:
            break
        arcos_caminho.append((atual, seguinte))
        atual = seguinte

    eh_valido, _, custo_verificado = verificar_caminho_valido(
        arcos_caminho, n, s=s, t=t, matriz_custos=c
    )

    gap_relativo = getattr(res, "mip_gap", 0.0)
    if gap_relativo is None or np.isnan(gap_relativo):
        gap_relativo = 0.0

    status_str = "Otimo" if res.status == 0 else "Limite de Tempo"

    return {
        "sucesso": True,
        "status": status_str,
        "valor_objetivo": float(res.fun),
        "gap_percentual": float(gap_relativo) * 100.0,
        "tempo_segundos": tempo,
        "arcos": arcos_caminho,
        "caminho_valido": eh_valido,
        "custo_verificado": custo_verificado,
    }


def main():
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
        print(f"Nenhuma instancia encontrada em: {pasta_instancias}")
        return

    resultados_gerais = []

    print("-" * 100)
    print(
        f"{'Instancia':<12} | {'n':>3} | {'Relax. Linear':>14} | {'Tempo LP':>9} | "
        f"{'Sol. Inteira':>14} | {'Gap (%)':>8} | {'Tempo PLI':>9} | {'Valido'}"
    )
    print("-" * 100)

    for caminho in arquivos:
        nome_arq = os.path.basename(caminho)
        _, n, c = ler_instancia_atsp(caminho)
        s, t = 1, n

        res_lp = resolver_relaxacao_linear(n, c, s=s, t=t)
        res_pli = resolver_modelo_inteiro(n, c, s=s, t=t, limite_tempo_segundos=limite_tempo)

        lp_val = res_lp["valor_objetivo"]
        pli_val = res_pli["valor_objetivo"]
        gap_pli = res_pli["gap_percentual"]
        valido_str = "SIM" if res_pli["caminho_valido"] else "NAO"

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

    print("-" * 100)

    caminho_csv = os.path.join(pasta_base, "resultados_instancias.csv")
    with open(caminho_csv, "w", newline="", encoding="utf-8") as f_csv:
        writer = csv.DictWriter(f_csv, fieldnames=list(resultados_gerais[0].keys()))
        writer.writeheader()
        writer.writerows(resultados_gerais)


if __name__ == "__main__":
    main()

