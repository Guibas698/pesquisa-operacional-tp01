#ifndef READFILE_H_
#define READFILE_H_

#include <stdio.h>
#include <vector>

int readATSPFile(FILE* graphFile, int *n, std::vector<std::vector<double>>& c);

#endif
