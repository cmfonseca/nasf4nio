/* ils.c
 *
 * (C) 2019, 2020 Carlos M. Fonseca <cmfonsec@dei.uc.pt>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License, version 3, as
 * published by the Free Software Foundation.
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
#include <gsl/gsl_rng.h>
#include "ils.h"

extern gsl_rng *rng;    /* The single rng instance used by the whole code */

struct solverState {
    struct problem *p;
    struct solution *s0, *s1;
    struct move *v;
    int kick_steps;
};

struct solverState *newSolver(struct problem *p) {
    struct solverState *ss;
    if (getNumObjectives(p) != 1) {
        fprintf(stderr, "ils error: problem must have a single objective.");
        return NULL;
    }
    ss = (struct solverState *) malloc(sizeof (struct solverState));
    ss->p = p;
    ss->s0 = allocSolution(p);
    ss->s1 = allocSolution(p);
    ss->v = allocMove(p);
    randomSolution(ss->s0);
    copySolution(ss->s1, ss->s0);
    /* default parameters */
    ss->kick_steps = 3;
    return ss;
}

void freeSolver(struct solverState *ss) {
    freeSolution(ss->s0);
    freeSolution(ss->s1);
    freeMove(ss->v);
    free(ss);
}

struct solverState *nextSolverState(struct solverState *ss) {
    struct solution *tmp;
    double o0, o1;
    if (randomMoveWOR(ss->v, ss->s1)) {
        if ((*getObjectiveIncrement(&o1, ss->v, ss->s1)) < 0)
            applyMove(ss->s1, ss->v);
    } else {
        if (*getObjectiveVector(&o1, ss->s1) <= *getObjectiveVector(&o0, ss->s0)) {
            tmp = ss->s0;
            ss->s0 = ss->s1;
            ss->s1 = tmp;
        }
        copySolution(ss->s1, ss->s0);
        /* printf("Kick!\n"); */
#if 1
        randomNeighbour(ss->s1, ss->kick_steps);
#else
        for (int i = 0; i < ss->kick_steps; i++)
            applyMove(ss->s1, randomMove(ss->v, ss->s1));
#endif
    }
    return ss;
}

struct solution *getSolverSolution(struct solverState *ss) {
    double o0, o1;
    if (*getObjectiveVector(&o1, ss->s1) < *getObjectiveVector(&o0, ss->s0))
        return ss->s1;
    else
        return ss->s0;
}

void printSolverState(struct solverState *ss) {
    printSolution(ss->s0);
    printSolution(ss->s1);
}

