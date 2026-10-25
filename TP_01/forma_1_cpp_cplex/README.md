# Forma 1: Implementação em C++ com IBM ILOG CPLEX (`forma_1_cpp_cplex`)

Esta pasta contém a resolução das **Atividades 1, 2 e 3** utilizando **C++ e a biblioteca IBM ILOG CPLEX Concert Technology**, seguindo exatamente o padrão do código-base fornecido pelo Prof. Paulo Henrique Macêdo na pasta `Implementacao/`.

## Por que o código original da pasta `Implementacao/` foi adaptado?
1. O arquivo original `Implementacao/readFile.c` lia grafos não-direcionados no formato DIMACS (`.col`, com linhas `p edge n m` e `e i j`), que possuía inclusive o `enum GT { graph, complement }` (usado em problemas de clique/coloração no grafo complemento).
2. Já as 13 instâncias disponibilizadas para este trabalho na pasta `Implementacao/instancias/Completos e Assimetricos nos Arcos/` são do formato **TSPLIB ATSP (`FULL_MATRIX`)** (grafos direcionados completos e assimétricos nos arcos, com matriz explícita de custos $c_{ij}$).
3. Portanto, adaptamos `readFile.h` e `readFile.c` para ler diretamente a matriz $n \times n$ das instâncias `.atsp` e implementamos o modelo completo em `PL.cpp`.

## Estrutura dos Arquivos
- `PL.cpp`:
  - **Atividade 1 (`verificarCaminhoValido`)**: Recebe uma sequência de arcos, verifica se tem $n-1$ arcos, se inicia em $s=1$ e termina em $t=n$, se respeita os graus de entrada/saída de todos os vértices e se percorre todos os $n$ vértices sem formar ciclos/subrotas.
  - **Atividade 2 (`resolverModeloMTZ(..., relaxacaoLinear = true)`)**: Resolve a **Relaxação Linear** do modelo com restrições de Miller-Tucker-Zemlin (Equações 1 a 11 do PDF), com variáveis contínuas `ILOFLOAT`.
  - **Atividade 3 (`resolverModeloMTZ(..., relaxacaoLinear = false)`)**: Adiciona a restrição de integralidade (`ILOBOOL` para $x_{ij}$ e `ILOINT` para $u_i$), configura o limite de 3 minutos (`cplex.setParam(IloCplex::Param::TimeLimit, 180.0)`), extrai o caminho encontrado e o valida chamando `verificarCaminhoValido`.
- `readFile.h` / `readFile.c`: Leitor do formato TSPLIB `.atsp` (`FULL_MATRIX`).
- `Makefile`: Configurado para compilação com `g++` e `CPLEX` no Linux/WSL.
- `executar_todas.sh`: Executa todas as 13 instâncias em sequência.
