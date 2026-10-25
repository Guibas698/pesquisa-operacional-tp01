# Trabalho Prático 01 - Pesquisa Operacional (QXD0181)

Implementação do modelo de Programação Linear Inteira (PLI) com restrições de Miller-Tucker-Zemlin (MTZ) para o Problema do Caminho Hamiltoniano Mínimo (PCHM) com vértices de origem ($s = 1$) e destino ($t = n$).

## Estrutura das Pastas

- `Implementacao/`: Código base disponibilizado na disciplina e pasta `instancias/Completos e Assimetricos nos Arcos/` contendo as 13 instâncias TSPLIB (`.atsp`).
- `TP_01/forma_1_cpp_cplex/`: Implementação em C++ seguindo a estrutura do código base (`PL.cpp`, `readFile.c`, `readFile.h`, `Makefile`).
- `TP_01/forma_2_python_mip/`: Implementação em Python (`pchm_mtz.py`) utilizando o solver HiGHS via `scipy.optimize` e arquivo `resultados_instancias.csv` com os dados coletados.
- `TP_01/RELATORIO_E_RESULTADOS.md`: Relatório referente à Atividade 4 com a comparação e análise dos resultados.

## Como Executar

### Versão C++ (`TP_01/forma_1_cpp_cplex`)
Compilação com `g++`:
```bash
g++ -O2 -static -I"TP_01/forma_1_cpp_cplex/include" "TP_01/forma_1_cpp_cplex/PL.cpp" -o "TP_01/forma_1_cpp_cplex/PL.exe"
```
Execução em uma instância:
```bash
./TP_01/forma_1_cpp_cplex/PL.exe "Implementacao/instancias/Completos e Assimetricos nos Arcos/br17.atsp"
```

### Versão Python (`TP_01/forma_2_python_mip`)
Execução em todas as instâncias (com limite padrão de 180s por instância):
```bash
python "TP_01/forma_2_python_mip/pchm_mtz.py"
```
Execução em uma instância específica:
```bash
python "TP_01/forma_2_python_mip/pchm_mtz.py" "Implementacao/instancias/Completos e Assimetricos nos Arcos/br17.atsp"
```
