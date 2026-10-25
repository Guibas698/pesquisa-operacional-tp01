# Trabalho PrÃ¡tico 01 - Pesquisa Operacional (QXD0181)

ImplementaÃ§Ã£o do modelo de ProgramaÃ§Ã£o Linear Inteira (PLI) com restriÃ§Ãµes de Miller-Tucker-Zemlin (MTZ) para o Problema do Caminho Hamiltoniano MÃ­nimo (PCHM) com vÃ©rtices de origem ($s = 1$) e destino ($t = n$).

## Estrutura das Pastas

- `Implementacao/`: CÃ³digo base disponibilizado na disciplina e pasta `instancias/Completos e Assimetricos nos Arcos/` contendo as 13 instÃ¢ncias TSPLIB (`.atsp`).
- `TP_01/forma_1_cpp_cplex/`: ImplementaÃ§Ã£o em C++ seguindo a estrutura do cÃ³digo base (`PL.cpp`, `readFile.c`, `readFile.h`, `Makefile`).
- `TP_01/forma_2_python_mip/`: ImplementaÃ§Ã£o em Python (`pchm_mtz.py`) utilizando o solver HiGHS via `scipy.optimize` e arquivo `resultados_instancias.csv` com os dados coletados.
- `TP_01/RELATORIO_E_RESULTADOS.md`: RelatÃ³rio referente Ã  Atividade 4 com a comparaÃ§Ã£o e anÃ¡lise dos resultados.

## Como Executar

### VersÃ£o C++ (`TP_01/forma_1_cpp_cplex`)
CompilaÃ§Ã£o com `g++`:
```bash
g++ -O2 -static -I"TP_01/forma_1_cpp_cplex/include" "TP_01/forma_1_cpp_cplex/PL.cpp" -o "TP_01/forma_1_cpp_cplex/PL.exe"
```
ExecuÃ§Ã£o em uma instÃ¢ncia:
```bash
./TP_01/forma_1_cpp_cplex/PL.exe "Implementacao/instancias/Completos e Assimetricos nos Arcos/br17.atsp"
```

### VersÃ£o Python (`TP_01/forma_2_python_mip`)
ExecuÃ§Ã£o em todas as instÃ¢ncias (com limite padrÃ£o de 180s por instÃ¢ncia):
```bash
python "TP_01/forma_2_python_mip/pchm_mtz.py"
```
ExecuÃ§Ã£o em uma instÃ¢ncia especÃ­fica:
```bash
python "TP_01/forma_2_python_mip/pchm_mtz.py" "Implementacao/instancias/Completos e Assimetricos nos Arcos/br17.atsp"
```

