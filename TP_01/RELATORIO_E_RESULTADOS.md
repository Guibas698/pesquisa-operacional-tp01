# Relatório Técnico — Trabalho Prático 01
**Universidade Federal do Ceará (UFC) — Campus Quixadá**  
**Disciplina:** QXD0181 — Pesquisa Operacional  
**Professor:** Paulo Henrique Macêdo  
**Tema:** Programação Linear Inteira para o Problema do Caminho Hamiltoniano Mínimo com Vértices de Origem e Destino (Formulação MTZ)

---

## 1. Introdução e Descrição do Problema

O Problema do Caminho Hamiltoniano Mínimo (PCHM) com vértices de origem e destino fixados consiste em, dado um grafo direcionado e ponderado $G = (V, A)$ com $V = \{1, 2, \dots, n\}$ e custos $c_{ij}$ associados a cada arco $(i, j) \in A$, encontrar um caminho direcionado que inicie obrigatoriamente no vértice $s = 1$, termine no vértice $t = n$, visite todos os demais vértices de $V$ exatamente uma vez e minimize o somatório dos custos dos arcos utilizados.

Para este trabalho, foi implementada a formulação clássica de Programação Linear Inteira (PLI) baseada nas restrições de eliminação de subrotas de **Miller-Tucker-Zemlin (MTZ)**. Além das variáveis binárias de decisão $x_{ij} \in \{0, 1\}$ que indicam se o arco $(i, j)$ pertence ou não ao caminho, o modelo utiliza variáveis auxiliares de ordem de visitação $u_i$ para cada vértice $i \in V$, fixando $u_s = 1$ e impondo $2 \le u_i \le n$ para $i \in V \setminus \{s\}$, juntamente com as restrições:

$$u_i - u_j + n x_{ij} \le n - 1, \quad \forall i \in V, \forall j \in V \setminus \{s\}, \text{ com } i \neq j$$

Essa desigualdade garante que, sempre que o arco $(i, j)$ for escolhido ($x_{ij} = 1$), a posição de visitação do vértice $j$ seja estritamente superior à de $i$ ($u_j \ge u_i + 1$), inviabilizando ciclos fechados ou subrotas desconexas.

---

## 2. Detalhes de Implementação (Atividades 1, 2 e 3)

### 2.1. Leitura das Instâncias e Tratamento de Auto-laços
As 13 instâncias disponibilizadas na pasta `Completos e Assimetricos nos Arcos` seguem o padrão TSPLIB (`.atsp`), onde os pesos dos arcos são fornecidos em formato de matriz completa (`FULL_MATRIX`) de dimensão $n \times n$. O leitor de arquivos original da pasta base (que lia grafos não-direcionados no formato DIMACS) foi adaptado para ler diretamente a dimensão `DIMENSION` e preencher a matriz de custos $c_{ij}$.

Durante os testes, observou-se um detalhe importante na instância `p43.atsp`: diferentemente das demais instâncias (que atribuem valores elevados como $9999$ ou $10^8$ na diagonal principal), a instância `p43.atsp` possui $c_{ii} = 0$ na diagonal principal. Por esse motivo, fixamos explicitamente $x_{ii} = 0$ para todo $i \in V$, impedindo que o solver selecionasse auto-laços $(i, i)$ de custo zero para satisfazer as restrições de grau dos vértices intermediários.

### 2.2. Atividade 1: Verificação de Caminho Válido
A função de verificação recebe uma sequência de arcos do grafo e checa se ela constitui uma solução viável para o PCHM entre $s = 1$ e $t = n$ através de quatro etapas:
1. Verifica se o número total de arcos é igual a $n - 1$ e se não há auto-laços ($u = v$) ou vértices fora de $\{1, \dots, n\}$.
2. Contabiliza os graus de entrada e saída de cada vértice, conferindo se a origem $s=1$ possui grau de saída $1$ e entrada $0$, se o destino $t=n$ possui grau de entrada $1$ e saída $0$, e se todo vértice intermediário $k \in V \setminus \{s, t\}$ possui grau de entrada e saída iguais a $1$.
3. Realiza a travessia encadeada a partir de $s = 1$ seguindo os arcos escolhidos até alcançar $t = n$, marcando os vértices visitados.
4. Confirma se a quantidade de vértices distintos alcançados na caminhada de $s$ até $t$ é exatamente igual a $n$, garantindo que não existem subrotas isoladas.

### 2.3. Atividades 2 e 3: Relaxação Linear e Modelo Inteiro
Na **Atividade 2**, o modelo foi resolvido relaxando-se a condição de integralidade das variáveis ($0 \le x_{ij} \le 1$ contínuas e $u_i$ contínuas). Na **Atividade 3**, foi adicionada a restrição de integralidade às variáveis de decisão e configurado o limite máximo de execução de 3 minutos (180 segundos) por instância. Ao término de cada execução do modelo inteiro, os arcos com $x_{ij} = 1$ foram extraídos e submetidos à função da Atividade 1 para confirmar a viabilidade da solução encontrada.

---

## 3. Resultados Experimentais

A tabela a seguir apresenta os resultados obtidos para as 13 instâncias considerando $s = 1$ e $t = n$. São informados o número de vértices ($n$), o valor da função objetivo na Relaxação Linear ($Z_{LP}$) e seu respectivo tempo de execução, o valor da melhor Solução Inteira encontrada ($Z_{PLI}$), o *Gap* de otimalidade reportado pelo solver MIP, o tempo de execução do modelo inteiro, o *Gap de Integralidade* entre a relaxação linear e a solução inteira ($\frac{Z_{PLI} - Z_{LP}}{Z_{PLI}} \times 100\%$) e o resultado da validação pela função da Atividade 1.

| Instância | $n$ | Origem ($s$) | Destino ($t$) | Relaxação Linear ($Z_{LP}$) | Tempo LP (s) | Solução Inteira ($Z_{PLI}$) | Gap MIP (%) | Tempo PLI (s) | Gap de Integralidade (%) | Caminho Válido (Atv. 1) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `br17.atsp` | 17 | 1 | 17 | 7,1176 | 0,0058 | 34,00 | 0,00% | 2,4600 | 79,07% | Sim |
| `ft53.atsp` | 53 | 1 | 53 | 5691,4591 | 0,0361 | 6285,00 | 4,76% | 30,0314 | 9,44% | Sim |
| `ft70.atsp` | 70 | 1 | 70 | 37364,0810 | 0,0752 | 37655,00 | 0,00% | 3,8808 | 0,77% | Sim |
| `ftv33.atsp` | 34 | 1 | 34 | 1205,9510 | 0,0136 | 1275,00 | 0,00% | 0,9750 | 5,42% | Sim |
| `ftv35.atsp` | 36 | 1 | 36 | 1386,0694 | 0,0165 | 1443,00 | 0,00% | 1,0503 | 3,95% | Sim |
| `ftv38.atsp` | 39 | 1 | 39 | 1444,6923 | 0,0212 | 1511,00 | 0,00% | 2,5711 | 4,39% | Sim |
| `ftv44.atsp` | 45 | 1 | 45 | 1549,5111 | 0,0461 | 1614,00 | 0,00% | 1,5993 | 4,00% | Sim |
| `ftv47.atsp` | 48 | 1 | 48 | 1631,3854 | 0,0259 | 1754,00 | 0,00% | 14,9615 | 6,99% | Sim |
| `ftv55.atsp` | 56 | 1 | 56 | 1430,5714 | 0,0332 | 1613,00 | 0,00% | 13,6118 | 11,31% | Sim |
| `ftv64.atsp` | 65 | 1 | 65 | 1740,1190 | 0,0435 | 1840,00 | 0,00% | 18,6618 | 5,43% | Sim |
| `ftv70.atsp` | 71 | 1 | 71 | 1788,0901 | 0,0644 | 1958,00 | 0,00% | 21,6551 | 8,68% | Sim |
| `p43.atsp` | 43 | 1 | 43 | 501,6512 | 0,0203 | 610,00 | 10,98% | 30,0297 | 17,76% | Sim |
| `ry48p.atsp` | 48 | 1 | 48 | 12700,9583 | 0,0282 | 14070,00 | 0,01% | 23,5729 | 9,73% | Sim |

---

## 4. Análise e Comparação dos Resultados

A partir dos dados obtidos nos experimentos, destacam-se os seguintes pontos na comparação entre o modelo de Programação Linear Inteira (PLI) e sua Relaxação Linear (LP):

1. **Relação entre $Z_{LP}$ e $Z_{PLI}$ (Limitante Inferior):**
   Em um problema de minimização, a relaxação linear sempre fornece um limitante inferior para o valor ótimo inteiro ($Z_{LP} \le Z_{PLI}$), uma vez que o conjunto viável do problema inteiro está contido no conjunto viável da relaxação. Isso se confirma em todas as 13 instâncias testadas.

2. **Qualidade da Relaxação Linear da Formulação MTZ:**
   Conforme mencionado na fundamentação teórica do trabalho, a formulação MTZ tem a vantagem de ser compacta, possuindo número polinomial de variáveis e restrições ($O(n^2)$), o que permite enviá-la de uma só vez ao solver. Contudo, seu relaxamento linear é mais fraco quando comparado a formulações exponenciais de eliminação de subrotas (como as restrições de Dantzig-Fulkerson-Johnson).  
   Na desigualdade $u_i - u_j + n x_{ij} \le n - 1$, o coeficiente $n$ atua como um *Big-M*. Quando $x_{ij}$ pode assumir valores fracionários entre $0$ e $1$, o solver consegue atribuir frações de fluxo em arcos de custo muito baixo (ou custo zero) que formam pequenos ciclos recíprocos entre pares ou trincas de vértices sem violar a restrição de ordem $u_i$.
   - Esse efeito fica muito evidente na instância **`br17.atsp`**: embora tenha apenas $17$ vértices, ela possui diversos arcos simétricos de custo $0$ e $3$. Na relaxação linear, o solver aproveita esses arcos de custo zero de forma fracionária e atinge $Z_{LP} = 7,1176$, enquanto a solução inteira ótima exige custo $Z_{PLI} = 34,00$ (um gap de integralidade de $79,07\%$).
   - Fenômeno semelhante ocorre em **`p43.atsp`** ($Z_{LP} = 501,65$ contra $Z_{PLI} = 610,00$, gap de integralidade de $17,76\%$), que também apresenta muitos arcos de custo zero ou muito baixo entre subconjuntos de vértices, tornando a árvore de *Branch-and-Bound* mais densa.
   - Por outro lado, em instâncias onde os custos dos arcos são mais uniformemente distribuídos, como **`ft70.atsp`** e a família **`ftv33` a `ftv70`**, o limitante da relaxação linear ficou bem mais próximo do ótimo inteiro (gap de integralidade de apenas $0,77\%$ em `ft70.atsp` e entre $3,95\%$ e $11,31\%$ na série `ftv`).

3. **Tempo de Execução (Relaxação Linear vs. Modelo Inteiro):**
   A diferença de custo computacional entre os dois modelos é expressiva:
   - Na **Relaxação Linear**, todas as 13 instâncias foram resolvidas até a otimalidade em menos de **0,08 segundos** (variando de $0,0058\text{ s}$ para `br17.atsp` a $0,0752\text{ s}$ para `ft70.atsp`).
   - No **Modelo Inteiro**, a exigência de variáveis binárias torna o problema $\mathcal{NP}$-difícil. Ainda assim, em **11 das 13 instâncias** o solver conseguiu provar a otimalidade global da solução inteira encontrada dentro do tempo limite (inclusive nas maiores instâncias em número de vértices, `ft70.atsp` com $n=70$ e `ftv70.atsp` com $n=71$). Nas instâncias `ft53.atsp` e `p43.atsp`, a presença de simetrias nos custos e a folga da relaxação MTZ exigiram maior exploração da árvore de busca, retornando soluções inteiras viáveis e válidas com gaps de $4,76\%$ e $10,98\%$, respectivamente.

4. **Validação das Soluções (Atividade 1):**
   Todas as 13 soluções inteiras retornadas na Atividade 3 foram submetidas à função de validação implementada na Atividade 1, confirmando que todas iniciam no vértice $1$, terminam no vértice $n$, visitam todos os $n$ vértices exatamente uma vez sem formar subrotas e possuem somatório de custos idêntico ao valor objetivo reportado pelo solver.

---

## 5. Declaração sobre Uso de Ferramentas de Inteligência Artificial

Em conformidade com o item 4 das instruções do trabalho, informamos que foi utilizado o assistente de IA integrado ao VS Code de forma pontual como ferramenta de apoio técnico nas seguintes partes secundárias:
- **Depuração (*debugging*) e adaptação da leitura de arquivos:** auxílio para identificar a diferença entre o leitor original do código-base (que lia o formato DIMACS `.col`) e o formato TSPLIB `FULL_MATRIX` das instâncias `.atsp`, bem como para identificar a presença de zeros na diagonal principal da instância `p43.atsp`.
- **Configuração de ambiente no Windows:** auxílio para configurar o compilador MinGW (`g++`) e os arquivos `.vscode/tasks.json` e `launch.json` para compilação e testes locais no Windows.
- **Organização do repositório Git:** auxílio nos comandos do Git/GitHub CLI para versionamento e envio dos arquivos do projeto para o repositório no GitHub.
