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

#define TIME_LIMIT_SECONDS 180.0

ILOSTLBEGIN
using namespace std;

// Questao 1: verifica se a sequencia de arcos forma um caminho hamiltoniano valido de s ate t
bool verificarCaminhoValido(
    const vector<pair<int, int>>& arcos,
    int n,
    int s,
    int t,
    const vector<vector<double>>& c,
    double* custoTotal = NULL
) {
    if (custoTotal != NULL) *custoTotal = 0.0;

    if (s < 1 || s > n || t < 1 || t > n || s == t) {
        return false;
    }

    if ((int)arcos.size() != n - 1) {
        return false;
    }

    vector<int> grauEntrada(n + 1, 0);
    vector<int> grauSaida(n + 1, 0);
    vector<int> proximo(n + 1, -1);
    double somaCustos = 0.0;

    for (size_t k = 0; k < arcos.size(); k++) {
        int u = arcos[k].first;
        int v = arcos[k].second;

        if (u < 1 || u > n || v < 1 || v > n || u == v) {
            return false;
        }

        grauSaida[u]++;
        grauEntrada[v]++;
        proximo[u] = v;
        somaCustos += c[u - 1][v - 1];
    }

    // Restricoes de grau na origem e no destino
    if (grauSaida[s] != 1 || grauEntrada[s] != 0) return false;
    if (grauEntrada[t] != 1 || grauSaida[t] != 0) return false;

    // Restricoes de grau nos vertices intermediarios
    for (int k = 1; k <= n; k++) {
        if (k == s || k == t) continue;
        if (grauEntrada[k] != 1 || grauSaida[k] != 1) {
            return false;
        }
    }

    // Percorre o caminho a partir de s para checar se todos os vertices sao visitados sem ciclos
    vector<bool> visitado(n + 1, false);
    int atual = s;
    int contVisitados = 0;

    while (atual != -1) {
        if (visitado[atual]) return false;
        visitado[atual] = true;
        contVisitados++;

        if (atual == t) break;
        atual = proximo[atual];
    }

    if (atual != t || contVisitados != n) {
        return false;
    }

    if (custoTotal != NULL) *custoTotal = somaCustos;
    return true;
}

// Questoes 2 e 3: monta e resolve o modelo MTZ (relaxacao linear ou modelo inteiro)
void resolverModeloMTZ(
    int n,
    int s,
    int t,
    const vector<vector<double>>& c,
    bool relaxacaoLinear
) {
    IloEnv env;
    try {
        IloModel model(env);
        IloCplex cplex(env);

        IloArray<IloNumVarArray> x(env, n);
        for (int i = 0; i < n; i++) {
            x[i] = IloNumVarArray(env, n);
            for (int j = 0; j < n; j++) {
                char nome[32];
                sprintf(nome, "x_%d_%d", i + 1, j + 1);
                if (i == j) {
                    x[i][j] = IloNumVar(env, 0.0, 0.0, ILOFLOAT, nome);
                } else if (relaxacaoLinear) {
                    x[i][j] = IloNumVar(env, 0.0, 1.0, ILOFLOAT, nome);
                } else {
                    x[i][j] = IloNumVar(env, 0.0, 1.0, ILOBOOL, nome);
                }
            }
        }

        // Variaveis u_i do MTZ
        IloNumVarArray u(env, n);
        for (int i = 0; i < n; i++) {
            char nome[32];
            sprintf(nome, "u_%d", i + 1);
            if (i == s) {
                u[i] = IloNumVar(env, 1.0, 1.0, relaxacaoLinear ? ILOFLOAT : ILOINT, nome);
            } else {
                u[i] = IloNumVar(env, 2.0, (double)n, relaxacaoLinear ? ILOFLOAT : ILOINT, nome);
            }
        }

        // Funcao objetivo (Eq. 2)
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

        // Eq. (3): saida de s igual a 1
        {
            IloExpr exp(env);
            for (int j = 0; j < n; j++) {
                if (j != s) exp += x[s][j];
            }
            restricoes.add(IloRange(env, 1.0, exp, 1.0));
            exp.end();
        }

        // Eq. (4): chegada em t igual a 1
        {
            IloExpr exp(env);
            for (int i = 0; i < n; i++) {
                if (i != t) exp += x[i][t];
            }
            restricoes.add(IloRange(env, 1.0, exp, 1.0));
            exp.end();
        }

        // Eq. (5): nao retorna para s
        {
            IloExpr exp(env);
            for (int i = 0; i < n; i++) {
                if (i != s) exp += x[i][s];
            }
            restricoes.add(IloRange(env, 0.0, exp, 0.0));
            exp.end();
        }

        // Eq. (6): nao sai de t
        {
            IloExpr exp(env);
            for (int j = 0; j < n; j++) {
                if (j != t) exp += x[t][j];
            }
            restricoes.add(IloRange(env, 0.0, exp, 0.0));
            exp.end();
        }

        // Eq. (7) e Eq. (8): conservacao de fluxo nos demais vertices
        for (int k = 0; k < n; k++) {
            if (k == s || k == t) continue;

            IloExpr expIn(env);
            for (int i = 0; i < n; i++) {
                if (i != k) expIn += x[i][k];
            }
            restricoes.add(IloRange(env, 1.0, expIn, 1.0));
            expIn.end();

            IloExpr expOut(env);
            for (int j = 0; j < n; j++) {
                if (j != k) expOut += x[k][j];
            }
            restricoes.add(IloRange(env, 1.0, expOut, 1.0));
            expOut.end();
        }

        // Eq. (10): restricoes MTZ de eliminacao de subrotas
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

        cplex.setOut(env.getNullStream());
        cplex.setWarning(env.getNullStream());

        if (!relaxacaoLinear) {
            cplex.setParam(IloCplex::Param::TimeLimit, TIME_LIMIT_SECONDS);
        }

        double inicio = cplex.getCplexTime();
        bool resolveu = cplex.solve();
        double tempoGasto = cplex.getCplexTime() - inicio;

        if (!resolveu) {
            printf("Inviavel ou sem solucao encontrada.\n");
        } else {
            double valorObj = cplex.getObjValue();
            if (relaxacaoLinear) {
                printf("[Relaxacao Linear]\n");
                printf("Valor Objetivo: %.4f | Tempo: %.4f s\n\n", valorObj, tempoGasto);
            } else {
                double gap = 100.0 * cplex.getMIPRelativeGap();
                printf("[Modelo Inteiro - PLI]\n");
                printf("Valor Objetivo: %.2f | Gap: %.2f%% | Tempo: %.4f s\n", valorObj, gap, tempoGasto);

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

                printf("Caminho: ");
                for (size_t k = 0; k < caminhoArcos.size(); k++) {
                    if (k == 0) printf("%d", caminhoArcos[k].first);
                    printf(" -> %d", caminhoArcos[k].second);
                }
                printf("\n");

                double custoVerificado = 0.0;
                bool valido = verificarCaminhoValido(caminhoArcos, n, s + 1, t + 1, c, &custoVerificado);
                printf("Caminho valido: %s (Custo = %.2f)\n\n", valido ? "SIM" : "NAO", custoVerificado);
            }
        }
    } catch (IloException& e) {
        cerr << "Erro CPLEX: " << e << endl;
    } catch (...) {
        cerr << "Erro na execucao." << endl;
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
        printf("Erro ao ler arquivo: %s\n", caminhoArquivo);
        return -1;
    }
    fclose(pFile);

    int s = 0;
    int t = n - 1;

    printf("------------------------------------------------------------\n");
    printf("Instancia: %s (n = %d, s = 1, t = %d)\n", caminhoArquivo, n, n);
    printf("------------------------------------------------------------\n");
    fflush(stdout);

    resolverModeloMTZ(n, s, t, c, true);
    fflush(stdout);

    resolverModeloMTZ(n, s, t, c, false);
    fflush(stdout);

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc >= 2) {
        if (resolverArquivoInstancia(argv[1]) != 0) {
            printf("Nao foi possivel abrir o arquivo: %s\n", argv[1]);
            return 1;
        }
        return 0;
    }

    const char* prefixos[] = {
        "../instancias/Completos e Assimetricos nos Arcos/",
        "../../instancias/Completos e Assimetricos nos Arcos/",
        "../../../Implementacao/instancias/Completos e Assimetricos nos Arcos/",
        "../../Implementacao/instancias/Completos e Assimetricos nos Arcos/",
        "../Implementacao/instancias/Completos e Assimetricos nos Arcos/",
        "Implementacao/instancias/Completos e Assimetricos nos Arcos/"
    };

    const char* nomesInstancias[] = {
        "br17.atsp", "ftv33.atsp", "ftv35.atsp", "ftv38.atsp",
        "p43.atsp", "ftv44.atsp", "ftv47.atsp", "ry48p.atsp",
        "ft53.atsp", "ftv55.atsp", "ftv64.atsp", "ft70.atsp", "ftv70.atsp"
    };
    const int totalInstancias = 13;

    string pastaEncontrada = "";
    for (int p = 0; p < 6; p++) {
        string teste = string(prefixos[p]) + "br17.atsp";
        FILE* f = fopen(teste.c_str(), "r");
        if (f != NULL) {
            fclose(f);
            pastaEncontrada = prefixos[p];
            break;
        }
    }

    if (pastaEncontrada.empty()) {
        printf("Uso: %s <arquivo_instancia.atsp>\n", argv[0]);
        return 1;
    }

    for (int i = 0; i < totalInstancias; i++) {
        string caminhoCompleto = pastaEncontrada + nomesInstancias[i];
        resolverArquivoInstancia(caminhoCompleto.c_str());
    }

    return 0;
}

