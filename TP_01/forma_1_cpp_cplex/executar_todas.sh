#!/bin/bash
# Script para compilar e executar o modelo C++/CPLEX em todas as 13 instancias ATSP

make compila

INST_DIR="../../Implementacao/instancias/Completos e Assimetricos nos Arcos"

for arq in "$INST_DIR"/*.atsp; do
    ./PL "$arq"
done
