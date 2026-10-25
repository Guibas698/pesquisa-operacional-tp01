#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"readFile.h"

int readFile(FILE * graphFile, const graphType t, int *n, int *m) {
	char               type  = ' ';
	char               linestr[100];
	char *             datastr;
	int                r, i, j;
	bool               namePrinted = false;

	*n = 0;
	*m = 0;

	//printf("Reading graph...\n");

	while (type != 'p') {
		type = fgetc(graphFile);
		if (type != EOF) {

			/* header */
			if (type == 'c') {
				datastr = fgets(linestr, 100, graphFile);
				if (datastr != NULL) {
                    if(!namePrinted) {
                        printf("%s", linestr);
                        namePrinted = true;
                    }
                } else return -1;
			}

			/* Vertices */
			if (type == 'p') {
				datastr = fgets(linestr, 100, graphFile);
				if (datastr == NULL)
					return -1;

				datastr = strtok(linestr," ");

				datastr = strtok(NULL," ");
				*n = atoi(datastr);

				datastr = strtok(NULL," ");
				*m = atoi(datastr);
				if (t == complement)
					*m = (((*n)*(*n) - *n) >> 1) - *m;
			}
		}
	}

	////
	// Graph variables
	////

	//printf("Graph with %d vertices and %d edges.\n", *n, *m);

	if (t == graph)
		setAllEdges(0);
	else {
		setAllEdges(0xFF);
		for (i = 0; i < *n; i++)
			removeEdge(i, i);
	}

	type = fgetc(graphFile);
	while (type != EOF) {
		/* Edges */
		if (type == 'e') {
			datastr = fgets(linestr, 100, graphFile);
			if (datastr == NULL)
				return -1;

			datastr = strtok(linestr," ");
			i = atol(datastr) - 1;

			datastr = strtok(NULL," ");
			j = atol(datastr) - 1;

			//printf("edge (%ld,%ld)\n",i,j);

			if (t == graph) {
				addEdge(i, j);
				addEdge(j, i);
			}
			else {
				removeEdge(i, j);
				removeEdge(j, i);
			}
		}
		else {
			datastr = fgets(linestr, 100, graphFile);
			if (datastr == NULL) return -1;
		}
		type = fgetc(graphFile);
	}

	double d = 200*( (*m)/ (double) (((*n)*((*n)-1))));
	printf("n=%d m=%d d=%.0lf\n",*n,*m,d);
	return 0;
}
