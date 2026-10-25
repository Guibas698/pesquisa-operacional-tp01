#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "readFile.h"

int readATSPFile(FILE* graphFile, int *n, std::vector<std::vector<double>>& c) {
    char word[256];
    *n = 0;

    // Procura por DIMENSION e EDGE_WEIGHT_SECTION no cabecalho TSPLIB
    while (fscanf(graphFile, "%255s", word) == 1) {
        if (strncmp(word, "DIMENSION", 9) == 0) {
            // Se for "DIMENSION:" le o numero logo em seguida; se for "DIMENSION" le ":" e depois o numero
            if (strchr(word, ':') == NULL) {
                fscanf(graphFile, "%255s", word); // consome ":"
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
