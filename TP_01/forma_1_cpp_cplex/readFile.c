#ifndef READFILE_C_INCLUDED_
#define READFILE_C_INCLUDED_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "readFile.h"

// Leitura das instancias no formato TSPLIB (.atsp)
inline int readATSPFile(FILE* graphFile, int *n, std::vector<std::vector<double>>& c) {
    char word[256];
    *n = 0;

    while (fscanf(graphFile, "%255s", word) == 1) {
        if (strncmp(word, "DIMENSION", 9) == 0) {
            if (strchr(word, ':') == NULL) {
                if (fscanf(graphFile, "%255s", word) != 1) return -1;
            }
            if (fscanf(graphFile, "%d", n) != 1) {
                return -1;
            }
        } else if (strcmp(word, "EDGE_WEIGHT_SECTION") == 0) {
            break;
        }
    }

    if (*n <= 0) {
        return -1;
    }

    c.assign(*n, std::vector<double>(*n, 0.0));

    for (int i = 0; i < *n; i++) {
        for (int j = 0; j < *n; j++) {
            double val;
            if (fscanf(graphFile, "%lf", &val) != 1) {
                return -1;
            }
            c[i][j] = val;
        }
    }

    return 0;
}

#endif
