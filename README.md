# QXD0181 - Pesquisa Operacional  Trabalho Prático 01

**Universidade Federal do Ceará (UFC)  Campus Quixadá**  
**Disciplina:** QXD0181 - Pesquisa Operacional  
**Professor:** Paulo Henrique Macêdo  
**Tema:** Programação Linear Inteira para o Problema do Caminho Hamiltoniano Mínimo (PCHM) com Vértices de Origem ($s = 1$) e Destino ($t = n$) utilizando as restrições de Miller-Tucker-Zemlin (MTZ).

---

## 📁 Estrutura do Repositório

```text
├── Implementacao/                     # Código-base original do professor e instâncias de teste
│   └── instancias/
│       └── Completos e Assimetricos nos Arcos/   # 13 instâncias TSPLIB (.atsp)
├── TP_01/
│   ├── TP 01.pdf                      # Especificação oficial do Trabalho Prático 01
│   ├── RELATORIO_E_RESULTADOS.md      # Relatório comparativo (Relaxação Linear vs. Modelo Inteiro)
│   ├── forma_1_cpp_cplex/             # Implementação em C++ com IBM ILOG CPLEX Concert Technology
│   │   ├── PL.cpp                     # Atividades 1, 2 e 3 em C++ / CPLEX
│   │   ├── readFile.c / readFile.h    # Leitor de instâncias TSPLIB ATSP (FULL_MATRIX)
│   │   ├── Makefile                   # Script de compilação C++ / CPLEX
│   │   └── README.md                  # Instruções da versão C++
│   └── forma_2_python_mip/            # Implementação em Python com Solver MIP/LP HiGHS (SciPy)
│       ├── pchm_mtz.py                # Atividades 1, 2 e 3 executáveis diretamente
│       ├── resultados_instancias.csv  # Resultados gerados para todas as 13 instâncias
│       └── README.md                  # Instruções da versão Python
└── .vscode/                           # Configurações prontas para execução com 1 clique no VS Code
```

---

## 🚀 Como Rodar e Testar no VS Code

### Opção 1: Usando o botão "Run and Debug" (`F5`) do VS Code
O projeto já inclui o arquivo `.vscode/launch.json` pré-configurado com o interpretador Python (`C:\Python313\python.exe`). Basta:
1. Abrir a pasta raiz `pesquisa operacional` no **VS Code**.
2. Pressionar `Ctrl+Shift+D` (aba **Run and Debug**) ou apertar **`F5`**.
3. Escolher uma das opções no menu superior:
   - **`TP 01: Testar apenas br17.atsp (Instantaneo)`**: Roda a Atividade 1, 2 e 3 na menor instância em ~2 segundos.
   - **`TP 01: Testar apenas ftv33.atsp (Instantaneo)`**: Roda em ~1 segundo.
   - **`TP 01: Rodar Todas as Instancias (Modo Rapido 30s)`**: Roda todas as 13 instâncias com limite de 30s por instância.
   - **`TP 01: Rodar Todas as Instancias (Limite Oficial 3 min / 180s)`**: Roda todas as 13 instâncias com o limite de 3 minutos especificado no PDF.

### Opção 2: Pelo Terminal Integrado do VS Code (`Ctrl + '`)

Para testar **uma única instância** rapidamente na mão:
```powershell
C:\Python313\python.exe "TP_01\forma_2_python_mip\pchm_mtz.py" "Implementacao\instancias\Completos e Assimetricos nos Arcos\br17.atsp"
```

Para rodar **todas as 13 instâncias** (com limite rápido de 30s por instância):
```powershell
C:\Python313\python.exe "TP_01\forma_2_python_mip\pchm_mtz.py" --limite=30
```

Para rodar **todas as 13 instâncias** com o limite oficial de **3 minutos (180s)** da Atividade 3:
```powershell
C:\Python313\python.exe "TP_01\forma_2_python_mip\pchm_mtz.py"
```

---

## 📊 Resumo dos Resultados Obtidos (Atividades 1, 2 e 3)

Considerando $s = 1$ (origem) e $t = n$ (destino) para cada instância:

| Instância | $n$ | Relaxação Linear ($Z_{LP}$) | Tempo LP (s) | Solução Inteira ($Z_{PLI}$) | Gap MIP (%) | Tempo PLI (s) | Validação (Atv. 1) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `br17.atsp` | 17 | 7,1176 | 0,0058 | **34,00** | 0,00% | 2,4600 | VÁLIDO |
| `ft53.atsp` | 53 | 5691,4591 | 0,0361 | **6285,00** | 4,76% | 30,0314 | VÁLIDO |
| `ft70.atsp` | 70 | 37364,0810 | 0,0752 | **37655,00** | 0,00% | 3,8808 | VÁLIDO |
| `ftv33.atsp` | 34 | 1205,9510 | 0,0136 | **1275,00** | 0,00% | 0,9750 | VÁLIDO |
| `ftv35.atsp` | 36 | 1386,0694 | 0,0165 | **1443,00** | 0,00% | 1,0503 | VÁLIDO |
| `ftv38.atsp` | 39 | 1444,6923 | 0,0212 | **1511,00** | 0,00% | 2,5711 | VÁLIDO |
| `ftv44.atsp` | 45 | 1549,5111 | 0,0461 | **1614,00** | 0,00% | 1,5993 | VÁLIDO |
| `ftv47.atsp` | 48 | 1631,3854 | 0,0259 | **1754,00** | 0,00% | 14,9615 | VÁLIDO |
| `ftv55.atsp` | 56 | 1430,5714 | 0,0332 | **1613,00** | 0,00% | 13,6118 | VÁLIDO |
| `ftv64.atsp` | 65 | 1740,1190 | 0,0435 | **1840,00** | 0,00% | 18,6618 | VÁLIDO |
| `ftv70.atsp` | 71 | 1788,0901 | 0,0644 | **1958,00** | 0,00% | 21,6551 | VÁLIDO |
| `p43.atsp` | 43 | 501,6512 | 0,0203 | **610,00** | 10,98% | 30,0297 | VÁLIDO |
| `ry48p.atsp` | 48 | 12700,9583 | 0,0282 | **14070,00** | 0,01% | 23,5729 | VÁLIDO |
