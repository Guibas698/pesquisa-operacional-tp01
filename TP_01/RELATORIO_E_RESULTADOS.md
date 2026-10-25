# RelatÃ³rio TÃ©cnico â€” Trabalho PrÃ¡tico 01
**Universidade Federal do CearÃ¡ (UFC) â€” Campus QuixadÃ¡**  
**Disciplina:** QXD0181 â€” Pesquisa Operacional  
**Professor:** Paulo Henrique MacÃªdo  
**Tema:** ProgramaÃ§Ã£o Linear Inteira para o Problema do Caminho Hamiltoniano MÃ­nimo com VÃ©rtices de Origem e Destino (FormulaÃ§Ã£o MTZ)

---

## 1. IntroduÃ§Ã£o e DescriÃ§Ã£o do Problema

O Problema do Caminho Hamiltoniano MÃ­nimo (PCHM) com vÃ©rtices de origem e destino fixados consiste em, dado um grafo direcionado e ponderado $G = (V, A)$ com $V = \{1, 2, \dots, n\}$ e custos $c_{ij}$ associados a cada arco $(i, j) \in A$, encontrar um caminho direcionado que inicie obrigatoriamente no vÃ©rtice $s = 1$, termine no vÃ©rtice $t = n$, visite todos os demais vÃ©rtices de $V$ exatamente uma vez e minimize o somatÃ³rio dos custos dos arcos utilizados.

Para este trabalho, foi implementada a formulaÃ§Ã£o clÃ¡ssica de ProgramaÃ§Ã£o Linear Inteira (PLI) baseada nas restriÃ§Ãµes de eliminaÃ§Ã£o de subrotas de **Miller-Tucker-Zemlin (MTZ)**. AlÃ©m das variÃ¡veis binÃ¡rias de decisÃ£o $x_{ij} \in \{0, 1\}$ que indicam se o arco $(i, j)$ pertence ou nÃ£o ao caminho, o modelo utiliza variÃ¡veis auxiliares de ordem de visitaÃ§Ã£o $u_i$ para cada vÃ©rtice $i \in V$, fixando $u_s = 1$ e impondo $2 \le u_i \le n$ para $i \in V \setminus \{s\}$, juntamente com as restriÃ§Ãµes:

$$u_i - u_j + n x_{ij} \le n - 1, \quad \forall i \in V, \forall j \in V \setminus \{s\}, \text{ com } i \neq j$$

Essa desigualdade garante que, sempre que o arco $(i, j)$ for escolhido ($x_{ij} = 1$), a posiÃ§Ã£o de visitaÃ§Ã£o do vÃ©rtice $j$ seja estritamente superior Ã  de $i$ ($u_j \ge u_i + 1$), inviabilizando ciclos fechados ou subrotas desconexas.

---

## 2. Detalhes de ImplementaÃ§Ã£o (Atividades 1, 2 e 3)

### 2.1. Leitura das InstÃ¢ncias e Tratamento de Auto-laÃ§os
As 13 instÃ¢ncias disponibilizadas na pasta `Completos e Assimetricos nos Arcos` seguem o padrÃ£o TSPLIB (`.atsp`), onde os pesos dos arcos sÃ£o fornecidos em formato de matriz completa (`FULL_MATRIX`) de dimensÃ£o $n \times n$. O leitor de arquivos original da pasta base (que lia grafos nÃ£o-direcionados no formato DIMACS) foi adaptado para ler diretamente a dimensÃ£o `DIMENSION` e preencher a matriz de custos $c_{ij}$.

Durante os testes, observou-se um detalhe importante na instÃ¢ncia `p43.atsp`: diferentemente das demais instÃ¢ncias (que atribuem valores elevados como $9999$ ou $10^8$ na diagonal principal), a instÃ¢ncia `p43.atsp` possui $c_{ii} = 0$ na diagonal principal. Por esse motivo, fixamos explicitamente $x_{ii} = 0$ para todo $i \in V$, impedindo que o solver selecionasse auto-laÃ§os $(i, i)$ de custo zero para satisfazer as restriÃ§Ãµes de grau dos vÃ©rtices intermediÃ¡rios.

### 2.2. Atividade 1: VerificaÃ§Ã£o de Caminho VÃ¡lido
A funÃ§Ã£o de verificaÃ§Ã£o recebe uma sequÃªncia de arcos do grafo e checa se ela constitui uma soluÃ§Ã£o viÃ¡vel para o PCHM entre $s = 1$ e $t = n$ atravÃ©s de quatro etapas:
1. Verifica se o nÃºmero total de arcos Ã© igual a $n - 1$ e se nÃ£o hÃ¡ auto-laÃ§os ($u = v$) ou vÃ©rtices fora de $\{1, \dots, n\}$.
2. Contabiliza os graus de entrada e saÃ­da de cada vÃ©rtice, conferindo se a origem $s=1$ possui grau de saÃ­da $1$ e entrada $0$, se o destino $t=n$ possui grau de entrada $1$ e saÃ­da $0$, e se todo vÃ©rtice intermediÃ¡rio $k \in V \setminus \{s, t\}$ possui grau de entrada e saÃ­da iguais a $1$.
3. Realiza a travessia encadeada a partir de $s = 1$ seguindo os arcos escolhidos atÃ© alcanÃ§ar $t = n$, marcando os vÃ©rtices visitados.
4. Confirma se a quantidade de vÃ©rtices distintos alcanÃ§ados na caminhada de $s$ atÃ© $t$ Ã© exatamente igual a $n$, garantindo que nÃ£o existem subrotas isoladas.

### 2.3. Atividades 2 e 3: RelaxaÃ§Ã£o Linear e Modelo Inteiro
Na **Atividade 2**, o modelo foi resolvido relaxando-se a condiÃ§Ã£o de integralidade das variÃ¡veis ($0 \le x_{ij} \le 1$ contÃ­nuas e $u_i$ contÃ­nuas). Na **Atividade 3**, foi adicionada a restriÃ§Ã£o de integralidade Ã s variÃ¡veis de decisÃ£o e configurado o limite mÃ¡ximo de execuÃ§Ã£o de 3 minutos (180 segundos) por instÃ¢ncia. Ao tÃ©rmino de cada execuÃ§Ã£o do modelo inteiro, os arcos com $x_{ij} = 1$ foram extraÃ­dos e submetidos Ã  funÃ§Ã£o da Atividade 1 para confirmar a viabilidade da soluÃ§Ã£o encontrada.

---

## 3. Resultados Experimentais

A tabela a seguir apresenta os resultados obtidos para as 13 instÃ¢ncias considerando $s = 1$ e $t = n$. SÃ£o informados o nÃºmero de vÃ©rtices ($n$), o valor da funÃ§Ã£o objetivo na RelaxaÃ§Ã£o Linear ($Z_{LP}$) e seu respectivo tempo de execuÃ§Ã£o, o valor da melhor SoluÃ§Ã£o Inteira encontrada ($Z_{PLI}$), o *Gap* de otimalidade reportado pelo solver MIP, o tempo de execuÃ§Ã£o do modelo inteiro, o *Gap de Integralidade* entre a relaxaÃ§Ã£o linear e a soluÃ§Ã£o inteira ($\frac{Z_{PLI} - Z_{LP}}{Z_{PLI}} \times 100\%$) e o resultado da validaÃ§Ã£o pela funÃ§Ã£o da Atividade 1.

| InstÃ¢ncia | $n$ | Origem ($s$) | Destino ($t$) | RelaxaÃ§Ã£o Linear ($Z_{LP}$) | Tempo LP (s) | SoluÃ§Ã£o Inteira ($Z_{PLI}$) | Gap MIP (%) | Tempo PLI (s) | Gap de Integralidade (%) | Caminho VÃ¡lido (Atv. 1) |
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

## 4. AnÃ¡lise e ComparaÃ§Ã£o dos Resultados

A partir dos dados obtidos nos experimentos, destacam-se os seguintes pontos na comparaÃ§Ã£o entre o modelo de ProgramaÃ§Ã£o Linear Inteira (PLI) e sua RelaxaÃ§Ã£o Linear (LP):

1. **RelaÃ§Ã£o entre $Z_{LP}$ e $Z_{PLI}$ (Limitante Inferior):**
   Em um problema de minimizaÃ§Ã£o, a relaxaÃ§Ã£o linear sempre fornece um limitante inferior para o valor Ã³timo inteiro ($Z_{LP} \le Z_{PLI}$), uma vez que o conjunto viÃ¡vel do problema inteiro estÃ¡ contido no conjunto viÃ¡vel da relaxaÃ§Ã£o. Isso se confirma em todas as 13 instÃ¢ncias testadas.

2. **Qualidade da RelaxaÃ§Ã£o Linear da FormulaÃ§Ã£o MTZ:**
   Conforme mencionado na fundamentaÃ§Ã£o teÃ³rica do trabalho, a formulaÃ§Ã£o MTZ tem a vantagem de ser compacta, possuindo nÃºmero polinomial de variÃ¡veis e restriÃ§Ãµes ($O(n^2)$), o que permite enviÃ¡-la de uma sÃ³ vez ao solver. Contudo, seu relaxamento linear Ã© mais fraco quando comparado a formulaÃ§Ãµes exponenciais de eliminaÃ§Ã£o de subrotas (como as restriÃ§Ãµes de Dantzig-Fulkerson-Johnson).  
   Na desigualdade $u_i - u_j + n x_{ij} \le n - 1$, o coeficiente $n$ atua como um *Big-M*. Quando $x_{ij}$ pode assumir valores fracionÃ¡rios entre $0$ e $1$, o solver consegue atribuir fraÃ§Ãµes de fluxo em arcos de custo muito baixo (ou custo zero) que formam pequenos ciclos recÃ­procos entre pares ou trincas de vÃ©rtices sem violar a restriÃ§Ã£o de ordem $u_i$.
   - Esse efeito fica muito evidente na instÃ¢ncia **`br17.atsp`**: embora tenha apenas $17$ vÃ©rtices, ela possui diversos arcos simÃ©tricos de custo $0$ e $3$. Na relaxaÃ§Ã£o linear, o solver aproveita esses arcos de custo zero de forma fracionÃ¡ria e atinge $Z_{LP} = 7,1176$, enquanto a soluÃ§Ã£o inteira Ã³tima exige custo $Z_{PLI} = 34,00$ (um gap de integralidade de $79,07\%$).
   - FenÃ´meno semelhante ocorre em **`p43.atsp`** ($Z_{LP} = 501,65$ contra $Z_{PLI} = 610,00$, gap de integralidade de $17,76\%$), que tambÃ©m apresenta muitos arcos de custo zero ou muito baixo entre subconjuntos de vÃ©rtices, tornando a Ã¡rvore de *Branch-and-Bound* mais densa.
   - Por outro lado, em instÃ¢ncias onde os custos dos arcos sÃ£o mais uniformemente distribuÃ­dos, como **`ft70.atsp`** e a famÃ­lia **`ftv33` a `ftv70`**, o limitante da relaxaÃ§Ã£o linear ficou bem mais prÃ³ximo do Ã³timo inteiro (gap de integralidade de apenas $0,77\%$ em `ft70.atsp` e entre $3,95\%$ e $11,31\%$ na sÃ©rie `ftv`).

3. **Tempo de ExecuÃ§Ã£o (RelaxaÃ§Ã£o Linear vs. Modelo Inteiro):**
   A diferenÃ§a de custo computacional entre os dois modelos Ã© expressiva:
   - Na **RelaxaÃ§Ã£o Linear**, todas as 13 instÃ¢ncias foram resolvidas atÃ© a otimalidade em menos de **0,08 segundos** (variando de $0,0058\text{ s}$ para `br17.atsp` a $0,0752\text{ s}$ para `ft70.atsp`).
   - No **Modelo Inteiro**, a exigÃªncia de variÃ¡veis binÃ¡rias torna o problema $\mathcal{NP}$-difÃ­cil. Ainda assim, em **11 das 13 instÃ¢ncias** o solver conseguiu provar a otimalidade global da soluÃ§Ã£o inteira encontrada dentro do tempo limite (inclusive nas maiores instÃ¢ncias em nÃºmero de vÃ©rtices, `ft70.atsp` com $n=70$ e `ftv70.atsp` com $n=71$). Nas instÃ¢ncias `ft53.atsp` e `p43.atsp`, a presenÃ§a de simetrias nos custos e a folga da relaxaÃ§Ã£o MTZ exigiram maior exploraÃ§Ã£o da Ã¡rvore de busca, retornando soluÃ§Ãµes inteiras viÃ¡veis e vÃ¡lidas com gaps de $4,76\%$ e $10,98\%$, respectivamente.

4. **ValidaÃ§Ã£o das SoluÃ§Ãµes (Atividade 1):**
   Todas as 13 soluÃ§Ãµes inteiras retornadas na Atividade 3 foram submetidas Ã  funÃ§Ã£o de validaÃ§Ã£o implementada na Atividade 1, confirmando que todas iniciam no vÃ©rtice $1$, terminam no vÃ©rtice $n$, visitam todos os $n$ vÃ©rtices exatamente uma vez sem formar subrotas e possuem somatÃ³rio de custos idÃªntico ao valor objetivo reportado pelo solver.

---

## 5. DeclaraÃ§Ã£o sobre Uso de Ferramentas de InteligÃªncia Artificial

Em conformidade com o item 4 das instruÃ§Ãµes do trabalho, informamos que foi utilizado o assistente de IA integrado ao VS Code de forma pontual como ferramenta de apoio tÃ©cnico nas seguintes partes secundÃ¡rias:
- **DepuraÃ§Ã£o (*debugging*) e adaptaÃ§Ã£o da leitura de arquivos:** auxÃ­lio para identificar a diferenÃ§a entre o leitor original do cÃ³digo-base (que lia o formato DIMACS `.col`) e o formato TSPLIB `FULL_MATRIX` das instÃ¢ncias `.atsp`, bem como para identificar a presenÃ§a de zeros na diagonal principal da instÃ¢ncia `p43.atsp`.
- **ConfiguraÃ§Ã£o de ambiente no Windows:** auxÃ­lio para configurar o compilador MinGW (`g++`) e os arquivos `.vscode/tasks.json` e `launch.json` para compilaÃ§Ã£o e testes locais no Windows.
- **OrganizaÃ§Ã£o do repositÃ³rio Git:** auxÃ­lio nos comandos do Git/GitHub CLI para versionamento e envio dos arquivos do projeto para o repositÃ³rio no GitHub.

