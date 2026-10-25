#ifndef READFILE_H_
#define READFILE_H_

#include <stdio.h>
#include <vector>

// Le um arquivo no formato TSPLIB ATSP (FULL_MATRIX)
// Preenche o numero de vertices n e a matriz de custos c (indexada de 0 a n-1)
int readATSPFile(FILE* graphFile, int *n, std::vector<std::vector<double>>& c);

#endif /*READFILE_H_*/
