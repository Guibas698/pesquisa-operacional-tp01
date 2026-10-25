// ============================================================================
// UFC - Campus Quixada
// QXD0181 - Pesquisa Operacional - Trabalho Pratico 01
// Prof. Paulo Henrique Macedo
//
// Programacao Linear Inteira para o Problema do Caminho Hamiltoniano Minimo
// com Vertices de Origem (s = 1) e Destino (t = n) - Formulacao MTZ
//
// Implementacao Forma 1: C++ com IBM ILOG CPLEX Concert Technology
// (Compativel tanto com Linux/CPLEX quanto com Windows/GCC no VS Code)
// ============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <cmath>

#if defined(__has_include)
  #if __has_include(<ilcplex/ilocplex.h>)
    #include <ilcplex/ilocplex.h>
  #else
    #include "include/ilcplex/ilocplex.h"
  #endif
#else
  #include <ilcplex/ilocplex.h>
#endif

#include "readFile.h"
#include "readFile.c"

#define EPS 1e-5
#define TIME_LIMIT_SECONDS 180.0 // 3 minutos (Atividade 3)

ILOSTLBEGIN
using namespace std;

// ============================================================================
// ATIVIDADE 1:
// Funcao que recebe uma sequencia de arcos do grafo e verifica se ela
// corresponde a um caminho valido, sendo uma solucao viavel para o PCHM
// com Vertices de Origem (s) e Destino (t).
//
// Nota: Os vertices nos arcos de entrada sao indexados de 1 a n (conforme PDF).
// ============================================================================
bool verificarCaminhoValido(
    const vector<pair<int, int>>& arcos,
    int n,
    int s,
    int t,
    const vector<vector<double>>& c,
    double* custoTotal = NULL,
    bool exibirDetalhes = false
) {
    if (custoTotal != NULL) *custoTotal = 0.0;

    // 1. Verifica se s e t sao validos e distintos
    if (s < 1 || s > n || t < 1 || t > n || s == t) {
        if (exibirDetalhes) printf("[Invalido] Vertices de origem (s=%d) ou destino (t=%d) invalidos.\n", s, t);
        return false;
    }

    // 2. Um caminho hamiltoniano em G=(V,A) com |V|=n deve ter exatamente n-1 arcos
    if ((int)arcos.size() != n - 1) {
        if (exibirDetalhes) printf("[Invalido] Numero de arcos (%d) diferente de n - 1 (%d).\n", (int)arcos.size(), n - 1);
        return false;
    }

    vector<int> grauEntrada(n + 1, 0);
    vector<int> grauSaida(n + 1, 0);
    vector<int> proximo(n + 1, -1);
    double somaCustos = 0.0;

    // 3. Verifica validade de cada arco e contabiliza graus de entrada e saida
    for (size_t k = 0; k < arcos.size(); k++) {
        int u = arcos[k].first;
        int v = arcos[k].second;

        // Vertices fora do intervalo {1, ..., n} ou auto-laco (u == v)
        if (u < 1 || u > n || v < 1 || v > n || u == v) {
            if (exibirDetalhes) printf("[Invalido] Arco (%d, %d) invalido ou auto-laco.\n", u, v);
            return false;
        }

        grauSaida[u]++;
        grauEntrada[v]++;
        proximo[u] = v;
        somaCustos += c[u - 1][v - 1];
    }

    // 4. Verifica as restricoes de grau (Equacoes 3 a 8 do PDF):
    // - Origem s: grau de saida = 1 (Eq. 3), grau de entrada = 0 (Eq. 5)
    if (grauSaida[s] != 1 || grauEntrada[s] != 0) {
        if (exibirDetalhes) printf("[Invalido] Violacao de grau na origem s=%d (saida=%d, entrada=%d).\n", s, grauSaida[s], grauEntrada[s]);
        return false;
    }

    // - Destino t: grau de entrada = 1 (Eq. 4), grau de saida = 0 (Eq. 6)
    if (grauEntrada[t] != 1 || grauSaida[t] != 0) {
        if (exibirDetalhes) printf("[Invalido] Violacao de grau no destino t=%d (entrada=%d, saida=%d).\n", t, grauEntrada[t], grauSaida[t]);
        return false;
    }

    // - Demais vertices k in V \ {s, t}: grau de entrada = 1 (Eq. 7) e saida = 1 (Eq. 8)
    for (int k = 1; k <= n; k++) {
        if (k == s || k == t) continue;
        if (grauEntrada[k] != 1 || grauSaida[k] != 1) {
            if (exibirDetalhes) printf("[Invalido] Violacao de grau no vertice intermediario k=%d (entrada=%d, saida=%d).\n", k, grauEntrada[k], grauSaida[k]);
            return false;
        }
    }

    // 5. Percorre o caminho a partir de s para garantir conectividade, ordem
    //    de visitacao unica e ausencia de subrotas disjuntas (ciclos isolados)
    vector<bool> visitado(n + 1, false);
    int atual = s;
    int contVisitados = 0;

    while (atual != -1) {
        if (visitado[atual]) {
            if (exibirDetalhes) printf("[Invalido] Ciclo detectado no vertice %d.\n", atual);
            return false;
        }
        visitado[atual] = true;
        contVisitados++;

        if (atual == t) break;
        atual = proximo[atual];
    }

    // Verifica se terminou em t e se visitou exatamente todos os n vertices
    if (atual != t || contVisitados != n) {
        if (exibirDetalhes) printf("[Invalido] O caminho de s=%d a t=%d visitou apenas %d de %d vertices (ha subrotas desconexas).\n", s, t, contVisitados, n);
        return false;
    }

    // Se a sequencia de arcos foi passada em ordem, tambem verifica continuidade direta
    bool sequenciaContigua = true;
    if (arcos[0].first != s || arcos.back().second != t) {
        sequenciaContigua = false;
    } else {
        for (size_t k = 0; k + 1 < arcos.size(); k++) {
            if (arcos[k].second != arcos[k + 1].first) {
                sequenciaContigua = false;
                break;
            }
        }
    }
    if (!sequenciaContigua && exibirDetalhes) {
        printf("[Aviso] Os arcos formam um Caminho Hamiltoniano valido, mas nao estavam ordenados sequencialmente no vetor.\n");
    }

    if (custoTotal != NULL) *custoTotal = somaCustos;
    return true;
}

// ============================================================================
// Construcao do Modelo MTZ no CPLEX (Atividades 2 e 3)
// Se relaxacaoLinear == true:
//   x_ij in [0, 1] continuas (IloNumVar) e u_i continuas (Atividade 2)
// Se relaxacaoLinear == false:
//   x_ij in {0, 1} binarias (IloBoolVar) e u_i inteiras (IloIntVar) (Atividade 3)
// ============================================================================
void resolverModeloMTZ(
    int n,
    int s, // indice 0-based (s = 0 corresponde ao vertice 1)
    int t, // indice 0-based (t = n-1 corresponde ao vertice n)
    const vector<vector<double>>& c,
    bool relaxacaoLinear
) {
    IloEnv env;
    try {
        IloModel model(env);
        IloCplex cplex(env);

        // Matriz de variaveis x_ij
        IloArray<IloNumVarArray> x(env, n);
        for (int i = 0; i < n; i++) {
            x[i] = IloNumVarArray(env, n);
            for (int j = 0; j < n; j++) {
                char nome[32];
                sprintf(nome, "x_%d_%d", i + 1, j + 1);
                if (i == j) {
                    // Auto-lacos nao fazem parte do caminho (x_ii = 0)
                    x[i][j] = IloNumVar(env, 0.0, 0.0, ILOFLOAT, nome);
                } else if (relaxacaoLinear) {
                    // Atividade 2: Relaxacao Linear (0 <= x_ij <= 1)
                    x[i][j] = IloNumVar(env, 0.0, 1.0, ILOFLOAT, nome);
                } else {
                    // Atividade 3: Variaveis Binarias x_ij in {0, 1} (Eq. 1)
                    x[i][j] = IloNumVar(env, 0.0, 1.0, ILOBOOL, nome);
                }
            }
        }

        // Variaveis de ordem de visitacao u_i (Secao 3 - MTZ)
        IloNumVarArray u(env, n);
        for (int i = 0; i < n; i++) {
            char nome[32];
            sprintf(nome, "u_%d", i + 1);
            if (i == s) {
                // Eq. (9): u_s = 1
                u[i] = IloNumVar(env, 1.0, 1.0, relaxacaoLinear ? ILOFLOAT : ILOINT, nome);
            } else {
                // Eq. (11): 2 <= u_i <= n, para todo i in V \ {s}
                u[i] = IloNumVar(env, 2.0, (double)n, relaxacaoLinear ? ILOFLOAT : ILOINT, nome);
            }
        }

        // Eq. (2): Funcao Objetivo: Minimizar Z = sum_{i in V} sum_{j in V} c_ij * x_ij
        IloExpr objExpr(env);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    objExpr += c[i][j] * x[i][j];
                }
            }
        }
        model.add(IloMinimize(env, objExpr));
        objExpr.end();

        IloRangeArray restricoes(env);

        // Eq. (3): sum_{j in V} x_{s,j} = 1 (O caminho deve sair de s)
        {
            IloExpr exp(env);
            for (int j = 0; j < n; j++) {
                if (j != s) exp += x[s][j];
            }
            restricoes.add(IloRange(env, 1.0, exp, 1.0, "saida_s"));
            exp.end();
        }

        // Eq. (4): sum_{i in V} x_{i,t} = 1 (O caminho deve chegar em t)
        {
            IloExpr exp(env);
            for (int i = 0; i < n; i++) {
                if (i != t) exp += x[i][t];
            }
            restricoes.add(IloRange(env, 1.0, exp, 1.0, "chegada_t"));
            exp.end();
        }

        // Eq. (5): sum_{i in V} x_{i,s} = 0 (O caminho nao pode retornar a s)
        {
            IloExpr exp(env);
            for (int i = 0; i < n; i++) {
                if (i != s) exp += x[i][s];
            }
            restricoes.add(IloRange(env, 0.0, exp, 0.0, "nao_retorna_s"));
            exp.end();
        }

        // Eq. (6): sum_{j in V} x_{t,j} = 0 (O caminho nao pode sair de t)
        {
            IloExpr exp(env);
            for (int j = 0; j < n; j++) {
                if (j != t) exp += x[t][j];
            }
            restricoes.add(IloRange(env, 0.0, exp, 0.0, "nao_sai_t"));
            exp.end();
        }

        // Eq. (7) e Eq. (8): Para todo k in V \ {s, t}
        for (int k = 0; k < n; k++) {
            if (k == s || k == t) continue;

            // Eq. (7): sum_{i in V} x_{i,k} = 1 (Chega em k apenas uma vez)
            IloExpr expIn(env);
            for (int i = 0; i < n; i++) {
                if (i != k) expIn += x[i][k];
            }
            restricoes.add(IloRange(env, 1.0, expIn, 1.0));
            expIn.end();

            // Eq. (8): sum_{j in V} x_{k,j} = 1 (Segue a partir de k apenas uma vez)
            IloExpr expOut(env);
            for (int j = 0; j < n; j++) {
                if (j != k) expOut += x[k][j];
            }
            restricoes.add(IloRange(env, 1.0, expOut, 1.0));
            expOut.end();
        }

        // Eq. (10): Restricoes MTZ: u_i - u_j + n * x_{i,j} <= n - 1
        // Para todo i in V, para todo j in V \ {s}, com i != j
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (j == s || i == j) continue;
                IloExpr expMTZ(env);
                expMTZ += u[i] - u[j] + n * x[i][j];
                restricoes.add(IloRange(env, -IloInfinity, expMTZ, n - 1));
                expMTZ.end();
            }
        }

        model.add(restricoes);
        cplex.extract(model);

        // Silencia saida padrao do CPLEX
        cplex.setOut(env.getNullStream());
        cplex.setWarning(env.getNullStream());

        // Limite de tempo de 3 minutos (180 segundos) para o modelo inteiro (Atividade 3)
        if (!relaxacaoLinear) {
            cplex.setParam(IloCplex::Param::TimeLimit, TIME_LIMIT_SECONDS);
        }

        double inicio = cplex.getCplexTime();
        bool resolveu = cplex.solve();
        double tempoGasto = cplex.getCplexTime() - inicio;

        if (!resolveu) {
            printf("Status: Inviavel ou sem solucao encontrada (Status = %d)\n", (int)cplex.getStatus());
        } else {
            double valorObj = cplex.getObjValue();
            if (relaxacaoLinear) {
                printf("=== ATIVIDADE 2: Relaxacao Linear (MTZ) ===\n");
                printf("Status          : %s\n", (cplex.getStatus() == IloAlgorithm::Optimal) ? "Otimo" : "Viavel");
                printf("Valor Objetivo  : %.4f\n", valorObj);
                printf("Tempo (s)       : %.4f s\n\n", tempoGasto);
            } else {
                double gap = 0.0;
                if (cplex.getStatus() == IloAlgorithm::Optimal || cplex.getStatus() == IloAlgorithm::Feasible) {
                    gap = 100.0 * cplex.getMIPRelativeGap();
                }
                printf("=== ATIVIDADE 3: Modelo Inteiro PLI (MTZ) ===\n");
                printf("Status          : %s\n", (cplex.getStatus() == IloAlgorithm::Optimal) ? "Otimo" : "Viavel (Limite de Tempo)");
                printf("Valor Objetivo  : %.2f\n", valorObj);
                printf("Gap Relativo    : %.2f%%\n", gap);
                printf("Tempo (s)       : %.4f s\n", tempoGasto);

                // Reconstroi a sequencia de arcos ordenada a partir de s
                vector<int> prox(n, -1);
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        if (i != j && cplex.getValue(x[i][j]) > 0.5) {
                            prox[i] = j;
                        }
                    }
                }

                vector<pair<int, int>> caminhoArcos;
                int atual = s;
                while (atual != -1 && atual != t && (int)caminhoArcos.size() < n) {
                    int seguinte = prox[atual];
                    if (seguinte == -1) break;
                    caminhoArcos.push_back(make_pair(atual + 1, seguinte + 1));
                    atual = seguinte;
                }

                // Exibe caminho e valida usando a funcao da Atividade 1
                printf("Caminho         : ");
                for (size_t k = 0; k < caminhoArcos.size(); k++) {
                    if (k == 0) printf("%d", caminhoArcos[k].first);
                    printf(" -> %d", caminhoArcos[k].second);
                }
                printf("\n");

                double custoVerificado = 0.0;
                bool valido = verificarCaminhoValido(caminhoArcos, n, s + 1, t + 1, c, &custoVerificado, true);
                printf("Validacao (Atv1): %s (Custo verificado = %.2f)\n\n",
                       valido ? "CAMINHO VALIDO (Solucao Viavel)" : "CAMINHO INVALIDO",
                       custoVerificado);
            }
        }
    } catch (IloException& e) {
        cerr << "Erro Concert/CPLEX: " << e << endl;
    } catch (...) {
        cerr << "Erro desconhecido." << endl;
    }
    env.end();
}

int resolverArquivoInstancia(const char* caminhoArquivo) {
    FILE* pFile = fopen(caminhoArquivo, "r");
    if (pFile == NULL) {
        return -1;
    }

    int n = 0;
    vector<vector<double>> c;
    if (readATSPFile(pFile, &n, c) != 0) {
        fclose(pFile);
        printf("Erro ao ler a instancia ATSP: %s\n", caminhoArquivo);
        return -1;
    }
    fclose(pFile);

    int s = 0;     // Vertice de origem: 1 (indice 0)
    int t = n - 1; // Vertice de destino: n (indice n-1)

    printf("============================================================\n");
    printf("Instancia: %s | n = %d | Origem s = 1 | Destino t = %d\n", caminhoArquivo, n, n);
    printf("============================================================\n");
    fflush(stdout);

    // Executa Atividade 2 (Relaxacao Linear)
    resolverModeloMTZ(n, s, t, c, true);
    fflush(stdout);

    // Executa Atividade 3 (Modelo Inteiro com limite de 3 minutos + Validacao da Atividade 1)
    resolverModeloMTZ(n, s, t, c, false);
    fflush(stdout);

    return 0;
}

int main(int argc, char* argv[]) {
    // Se um arquivo especifico foi passado via linha de comando, executa nele
    if (argc >= 2) {
        if (resolverArquivoInstancia(argv[1]) != 0) {
            printf("Nao foi possivel abrir o arquivo: %s\n", argv[1]);
            return 1;
        }
        return 0;
    }

    // Se executado diretamente pelo botao Play do VS Code (inclusive de dentro da pasta output/),
    // localiza automaticamente a pasta "Implementacao/instancias/Completos e Assimetricos nos Arcos"
    const char* prefixos[] = {
        "../../../Implementacao/instancias/Completos e Assimetricos nos Arcos/",
        "../../Implementacao/instancias/Completos e Assimetricos nos Arcos/",
        "../Implementacao/instancias/Completos e Assimetricos nos Arcos/",
        "Implementacao/instancias/Completos e Assimetricos nos Arcos/"
    };

    const char* nomesInstancias[] = {
        "br17.atsp",
        "ftv33.atsp",
        "ftv35.atsp",
        "ftv38.atsp",
        "p43.atsp",
        "ftv44.atsp",
        "ftv47.atsp",
        "ry48p.atsp",
        "ft53.atsp",
        "ftv55.atsp",
        "ftv64.atsp",
        "ft70.atsp",
        "ftv70.atsp"
    };
    const int totalInstancias = 13;

    string pastaEncontrada = "";
    for (int p = 0; p < 4; p++) {
        string teste = string(prefixos[p]) + "br17.atsp";
        FILE* f = fopen(teste.c_str(), "r");
        if (f != NULL) {
            fclose(f);
            pastaEncontrada = prefixos[p];
            break;
        }
    }

    if (pastaEncontrada.empty()) {
        printf("Nao foi possivel localizar a pasta de instancias automaticamente.\n");
        printf("Uso: %s <arquivo_instancia.atsp>\n", argv[0]);
        return 1;
    }

    for (int i = 0; i < totalInstancias; i++) {
        string caminhoCompleto = pastaEncontrada + nomesInstancias[i];
        resolverArquivoInstancia(caminhoCompleto.c_str());
    }

    return 0;
}
