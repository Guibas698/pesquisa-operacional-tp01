# Forma 2: Implementação em Python com Solver MIP/LP HiGHS (`forma_2_python_mip`)

Esta pasta contém a resolução completa e executável das **Atividades 1, 2 e 3** em **Python**, utilizando o solver de Programação Linear e Inteira Mista **HiGHS** (via `scipy.optimize.linprog` e `scipy.optimize.milp`).

## Como Executar

Para rodar em todas as 13 instâncias com o limite padrão de **3 minutos (180 segundos)** por instância (conforme exigido na Atividade 3):
```bash
python pchm_mtz.py
```

Para rodar em uma única instância específica:
```bash
python pchm_mtz.py "../../Implementacao/instancias/Completos e Assimetricos nos Arcos/br17.atsp"
```

Para rodar com um limite de tempo personalizado (ex.: 35 segundos por instância):
```bash
python pchm_mtz.py --limite=35
```

## O que o script faz:
1. **Atividade 1 (`verificar_caminho_valido`)**: Testa casos válidos e inválidos (com subrotas ou desconexos) e valida todas as soluções inteiras obtidas nas 13 instâncias.
2. **Atividade 2 (`resolver_relaxacao_linear`)**: Monta as Equações (2) a (11) do PDF e resolve a **Relaxação Linear** para cada instância `.atsp`.
3. **Atividade 3 (`resolver_modelo_inteiro`)**: Adiciona a restrição de integralidade ($x_{ij} \in \{0,1\}$ e $u_i \in \mathbb{Z}$) com limite de até 3 minutos, reconstrói a sequência de arcos do caminho de $s=1$ até $t=n$, verifica a solução com a função da Atividade 1 e salva a tabela completa em `resultados_instancias.csv`.
