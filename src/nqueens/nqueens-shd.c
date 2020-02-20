/* nqueens-shd.c
 *
 * (C) 2018 Carlos M. Fonseca <cmfonsec@dei.uc.pt>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include "nqueens.h"
#include "shd.h"

gsl_rng *rng;    /* The single rng instance used by the whole code */

int main(int argc, char **argv) {
    struct problem *p;
    struct solverState *ss;
    int board_size, max_iter, i;
    double cost, mincost;

    if (argc < 3) {
        fprintf(stderr, "Usage: %s <board size> <max iter>\n", argv[0]);
        return 0;
    }

    /* Set up random number generation */
    rng = gsl_rng_alloc(gsl_rng_mt19937); 
    gsl_rng_set(rng, time(0));

    /* Input arguments */
    board_size = atoi(argv[1]);
    max_iter = atoi(argv[2]);

    /* Problem and solver instantiation */
    p = newProblem(board_size);
    if (p != NULL) {
        ss = newSolver(p);

        /* Run */
        getObjectiveVector(&mincost, getSolverSolution(ss));
        printf("iter = 0, obj = %.0f\n", mincost);
        for (i = 0; mincost > 0 && i < max_iter; i++) {
            nextSolverState(ss);
            getObjectiveVector(&cost, getSolverSolution(ss));
            if (cost < mincost) {
                mincost = cost;
                printf("iter = %d, obj = %.0f\n", i+1, mincost);
            }
        }

        /* Report result */
        printSolution(getSolverSolution(ss));

        /* Clean up */
        freeSolver(ss);
        freeProblem(p);
    }
    gsl_rng_free(rng);
    return 0;
}

