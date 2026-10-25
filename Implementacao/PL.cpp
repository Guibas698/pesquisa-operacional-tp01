#include <stdio.h>
#include <time.h>
#include <iostream>
#include <limits.h>
#include <math.h>
#include <ilcplex/ilocplex.h>
#include <vector>
#include <list>
#include"readFile.h"

#define eps 0.000001

ILOSTLBEGIN
using namespace std;

int n,m;
IloEnv   env;
IloModel model(env);
IloNumVarArray var_x(env);
IloObjective obj = IloMaximize(env);
IloCplex cplex(env);

void createModel() {
    IloExpr exp(env);
    IloRangeArray con(env);

    try {
	double coef[n];

        //printf("Creating %d X variables...\n");
        for(int i=0;i<n;i++) {
                char nome[20];
                sprintf(nome,"x%d",i+1);
                var_x.add(IloNumVar(env,0,IloInfinity));
                var_x[i].setName(nome);
        }

        //printf("Creating Objective Function...\n\n");
        coef[0]=3;
        coef[1]=5;
        for(int i=0;i<n;i++) obj.setLinearCoef(var_x[i], coef[i]);

        //printf("Adding NO two classes constraints...\n\n");
	exp += var_x[0];
        con.add(IloRange(env, exp, 4));
        exp.clear();

        exp += 2*var_x[1];
        con.add(IloRange(env, exp, 12));
        exp.clear();

	exp += 3*var_x[0] + 2*var_x[1];
        con.add(IloRange(env, exp, 18));
        exp.clear();

        model.add(con);
        con.clear();
        model.add(obj);
        cplex.extract(model);
        cplex.exportModel("PL.lp");      // only for printing
    }
    catch (IloException& e) {
        cerr << "Concert exception caught: " << e << endl;
    }
    catch (...) {
        cerr << "Unknown exception caught" << endl;
    }
}

bool solveMaster(double *x, double *sol_val) {    // Return true if the master problem is feasible
    // Optimize the problem
    if ( !cplex.solve() && cplex.getStatus() != IloAlgorithm::Infeasible) {
        env.error() << "Failed to optimize LP" << endl;
        cout << cplex.getStatus() << endl;
        throw(-1);
    }

    // Get solution
    if(cplex.getStatus() != IloAlgorithm::Infeasible) {
        *sol_val = cplex.getObjValue();
        //printf("\nSolution value = %.2lf ",*sol_val);
        for(int i=0;i<n;i++)
            x[i] = cplex.getValue(var_x[i]);
        return true;
    }
    return false;
}

void Procedure() {
    double sol_val=-1;
    double *x=NULL;

    x=(double*) malloc((n)*sizeof(double));    
    for(int i=0;i<n;i++) memset(x,0.0,n*sizeof(double));

    if(!solveMaster(x,&sol_val)) printf("INFEASIBLE!\n");   // Return false if the master problem is infeasible

    printf("Objective Function Value = %.2lf\n",sol_val);
    for(int i=0;i<n;i++)
	printf("x[%d] = %2.lf\n", i+1, x[i]);

    free(x);
}


int main (int argc, char *argv[]) {
    if (argc < 2) {
        printf("Please, specify the file.\n");
        return 0;
    }
    FILE * pFile = fopen(argv[1], "r");
    if (pFile == NULL) {
        printf("Could not open %s.\n", argv[1]);
        return 0;
    }

    if (readFile(pFile, graph, &n, &m)) {
        fclose(pFile);
        printf("Could not read file %s.", argv[1]);
        return 0;
    }
    fclose(pFile);

    n=2;	// REMOVE THIS LINE TO START YOUR WORK!

    clock_t start, end;
    double elapsed;
    start = clock();
    // Solving the problem
    cplex.setOut(env.getNullStream());      // Do not print cplex informations
    cplex.setWarning(env.getNullStream());  // Do not print cplex warnings
    createModel();
    Procedure();

    // Free memory

    // Getting run time
    env.end();
    end = clock();
    elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Time: %.5g second(s).\n", elapsed);

    return 0;
}
