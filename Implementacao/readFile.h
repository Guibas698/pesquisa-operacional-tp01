#ifndef READFILE_H_
#define READFILE_H_

#include"graphMacros.h"
#include"bitMapGraph.h"
#include<stdio.h>

typedef enum GT { graph, complement } graphType;

int readFile(FILE* graphFile, const graphType t, int *n, int *m);

#endif /*READFILE_H_*/
