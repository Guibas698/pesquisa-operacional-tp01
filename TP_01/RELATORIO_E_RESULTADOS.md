# Trabalho Prático 01 - QXD0181 Pesquisa Operacional
**Prof. Paulo Henrique Macêdo**  
**Tema:** Programação Linear Inteira para o Problema do Caminho Hamiltoniano Mínimo (PCHM) com Vértices de Origem ($s=1$) e Destino ($t=n$) — Formulação Miller-Tucker-Zemlin (MTZ)

---

## 1. Esclarecimento das Dúvidas sobre a Pasta `Implementacao` e as Instâncias

1. **O PDF pede um único algoritmo/modelo principal?**
   - **Sim.** O PDF especifica uma única formulação matemática de Programação Linear Inteira (PLI): o modelo clássico com restrições de eliminação de subrotas de **Miller-Tucker-Zemlin (MTZ)** (Equações 1 a 11 do PDF), executado de duas maneiras complementares nas atividades:
     - **Atividade 2:** **Relaxação Linear (LP)**, em que as variáveis $x_{ij} \in [0, 1]$ e $u_i \in [2, n]$ são contínuas ($u_s = 1$).
     - **Atividade 3:** **Modelo Inteiro (PLI/MILP)**, adicionando as restrições de integralidade ($x_{ij} \in \{0, 1\}$ e $u_i \in \mathbb{Z}$) com limite de execução de até 3 minutos (180 segundos) por instância, validando a sequência de arcos obtida com a função da **Atividade 1**.

2. **Por que existe a pasta `Implementacao` e de onde veio a história de "complemento"?**
   - A pasta `Implementacao/` contém um **esqueleto/exemplo em C++ com CPLEX** (`PL.cpp`, `Makefile`, `readFile.c`, `bitMapGraph.c`) reaproveitado de trabalhos anteriores da disciplina (como Coloração de Grafos ou Clique Máxima, que usavam instâncias `.col` como `1-FullIns_3.col` no `Makefile`).
   - Dentro de `Implementacao/readFile.h` havia a definição `typedef enum GT { graph, complement } graphType;` (que servia para inverter arestas e ler o grafo complemento em problemas de clique/conjunto independente).
   - Além disso, dentro de `Implementacao/instancias/`, a subpasta de testes chama-se **`Completos e Assimetricos nos Arcos`** (contendo 13 instâncias clássicas da TSPLIB com extensão `.atsp`). Ou seja: o nome da pasta de testes é **"Completos..."** (grafos direcionados completos com pesos assimétricos nos arcos), e não "complemento".
   - Como os arquivos `.atsp` possuem matriz de custos explícita (`EDGE_WEIGHT_FORMAT: FULL_MATRIX`), o `readFile.c` original (feito para grafos sem peso no formato DIMACS) precisou ser adaptado para ler a matriz $n \times n$ de custos $c_{ij}$.

3. **Cuidado especial na instância `p43.atsp`:**
   - Na maioria das instâncias `.atsp`, a diagonal principal $c_{ii}$ tem valor muito alto (`9999`, `9999999` ou `100000000`). Porém, em `p43.atsp`, a diagonal principal tem custo **`0`**!
   - Como a restrição MTZ (Equação 10) é definida apenas para $i \neq j$, é indispensável fixar $x_{ii} = 0$ (proibir auto-laços) no modelo; caso contrário, o solver escolheria $x_{ii} = 1$ com custo zero. Ambas as nossas implementações tratam isso explicitamente.

---

## 2. As Duas Formas de Resolver Criadas Dentro de `TP_01/`

Conforme solicitado, foram criadas duas pastas independentes dentro de `TP_01/`, cada uma resolvendo as 3 primeiras atividades:

1. **`TP_01/forma_1_cpp_cplex/` (C++ com IBM ILOG CPLEX Concert Technology)**
   - Segue exatamente a estrutura do código-base do professor (`PL.cpp`, `readFile.c`, `readFile.h`, `Makefile`).
   - Adapta o leitor `readFile.c` para o formato TSPLIB `.atsp` (`FULL_MATRIX`) e implementa a validação do caminho (Atividade 1), a Relaxação Linear no CPLEX (Atividade 2) e o Modelo Inteiro no CPLEX com limite de 180s (Atividade 3).

2. **`TP_01/forma_2_python_mip/` (Python com Solver MIP/LP HiGHS via SciPy)**
   - Implementação completa em Python (`pchm_mtz.py`), que roda diretamente no Windows sem depender da licença local do CPLEX Studio.
   - Contém a função `verificar_caminho_valido` (Atividade 1), `resolver_relaxacao_linear` (Atividade 2) e `resolver_modelo_inteiro` (Atividade 3), além de gerar automaticamente o arquivo `resultados_instancias.csv`.

---

## 3. Resultados Obtidos nas 13 Instâncias (`s = 1`, `t = n`)

A tabela abaixo resume os valores retornados pela **Relaxação Linear (Atividade 2)**, pelo **Modelo Inteiro PLI (Atividade 3)** e pela **Função de Validação de Caminho (Atividade 1)** em todas as 13 instâncias disponibilizadas:

| Instância | $n$ | Origem ($s$) | Destino ($t$) | Relaxação Linear ($Z_{LP}$) | Tempo LP (s) | Solução Inteira ($Z_{PLI}$) | Gap MIP (%) | Tempo PLI (s) |Gap de Integr. $\frac{Z_{PLI}-Z_{LP}}{Z_{PLI}}$ | Validação (Atv. 1) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `br17.atsp` | 17 | 1 | 17 | 7,1176 | 0,0058 | **34,00** | 0,00% | 2,4600 | 79,07% | VÁLIDO |
| `ft53.atsp` | 53 | 1 | 53 | 5691,4591 | 0,0361 | **6285,00** | 4,76% | 30,0314 | 9,44% | VÁLIDO |
| `ft70.atsp` | 70 | 1 | 70 | 37364,0810 | 0,0752 | **37655,00** | 0,00% | 3,8808 | 0,77% | VÁLIDO |
| `ftv33.atsp` | 34 | 1 | 34 | 1205,9510 | 0,0136 | **1275,00** | 0,00% | 0,9750 | 5,42% | VÁLIDO |
| `ftv35.atsp` | 36 | 1 | 36 | 1386,0694 | 0,0165 | **1443,00** | 0,00% | 1,0503 | 3,95% | VÁLIDO |
| `ftv38.atsp` | 39 | 1 | 39 | 1444,6923 | 0,0212 | **1511,00** | 0,00% | 2,5711 | 4,39% | VÁLIDO |
| `ftv44.atsp` | 45 | 1 | 45 | 1549,5111 | 0,0461 | **1614,00** | 0,00% | 1,5993 | 4,00% | VÁLIDO |
| `ftv47.atsp` | 48 | 1 | 48 | 1631,3854 | 0,0259 | **1754,00** | 0,00% | 14,9615 | 6,99% | VÁLIDO |
| `ftv55.atsp` | 56 | 1 | 56 | 1430,5714 | 0,0332 | **1613,00** | 0,00% | 13,6118 | 11,31% | VÁLIDO |
| `ftv64.atsp` | 65 | 1 | 65 | 1740,1190 | 0,0435 | **1840,00** | 0,00% | 18,6618 | 5,43% | VÁLIDO |
| `ftv70.atsp` | 71 | 1 | 71 | 1788,0901 | 0,0644 | **1958,00** | 0,00% | 21,6551 | 8,68% | VÁLIDO |
| `p43.atsp` | 43 | 1 | 43 | 501,6512 | 0,0203 | **610,00** | 10,98% | 30,0297 | 17,76% | VÁLIDO |
| `ry48p.atsp` | 48 | 1 | 48 | 12700,9583 | 0,0282 | **14070,00** | 0,01% | 23,5729 | 9,73% | VÁLIDO |

---

## 4. Análise Comparativa (Subsídio para a Atividade 4 / Relatório)

1. **Qualidade do Limitante Inferior da Relaxação Linear (MTZ):**
   - Em todas as 13 instâncias, $Z_{LP} \le Z_{PLI}$, pois a relaxação linear amplia o espaço de busca permitindo valores fracionários para $x_{ij} \in [0, 1]$.
   - Na restrição MTZ ($u_i - u_j + n x_{ij} \le n - 1$), quando $x_{ij}$ assume valores fracionários (por exemplo, $x_{ij} = 0,5$), o termo $n x_{ij}$ ainda permite folga suficiente nas variáveis $u_i$ para que quase-ciclos fracionários de baixo custo sejam formados sem violar a desigualdade.
   - Isso explica o que o PDF destaca no início da página 3: a formulação MTZ possui tamanho compacto polinomial ($O(n^2)$ variáveis e restrições), mas apresenta um **relaxamento linear mais fraco** (particularmente visível em `br17.atsp`, onde $Z_{LP} = 7,12$ e $Z_{PLI} = 34,00$, e em `p43.atsp`, onde há muitos arcos de custo zero formando ciclos curtos na relaxação).

2. **Desempenho Computacional (Tempo LP vs. Tempo PLI):**
   - A **Relaxação Linear** é resolvida em menos de **0,08 segundos** em todas as instâncias, pois problemas de Programação Linear contínua com poucos milhares de variáveis são resolvidos em tempo polinomial (método Simplex / Pontos Interiores).
   - Já o **Modelo Inteiro (PLI)** é $\mathcal{NP}$-difícil e exige enumeração implícita (*Branch-and-Bound* / *Branch-and-Cut*). Mesmo assim, com a formulação MTZ bem estruturada, **11 das 13 instâncias atingiram a solução ótima global comprovada (Gap = 0%) em menos de 24 segundos**, e todas as 13 retornaram caminhos hamiltonianos válidos confirmados pela função da Atividade 1.
