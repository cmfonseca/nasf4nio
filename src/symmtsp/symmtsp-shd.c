
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include "symmtsp.h"
#include "shd.h"

gsl_rng *rng;    /* The single rng instance used by the whole code */

int main(int argc, char **argv) {
    struct problem *p;
    struct solverState *ss;
    int max_iter, i;
    double cost, mincost;
    char *filename;

    if (argc < 3) {
        fprintf(stderr, "Usage: %s <problem file> <max iter>\n", argv[0]);
        return 0;
    }

    /* Set up random number generation */
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, time(0));

    filename = argv[1];
    max_iter = atoi(argv[2]);

    p = newProblem(filename);
    if (p != NULL) {
        ss = newSolver(p);

        /* Run */
        getObjectiveVector(&mincost, getSolverSolution(ss));
        setObjValue(mincost, getSolverSolution(ss));
        printf("iter = 0, obj = %.0f\n", mincost);
        for (i = 0; mincost > 0 && i < max_iter; i++) {
            nextSolverState(ss);
            getObjectiveVector(&cost, getSolverSolution(ss));

            if (cost < mincost) {
            	setObjValue(cost,getSolverSolution(ss));
                mincost = cost;
                printf("iter = %d, obj = %.2lf\n", i+1, mincost);
            }
        }

        /* Report result */
        setObjValue(mincost, getSolverSolution(ss));
        printSolution(getSolverSolution(ss));
        /* Clean up */
        freeSolver(ss);
        freeProblem(p);
    }
    gsl_rng_free(rng);
    return 0;
}

