/* shd.c
 *
 * (C) 2018, 2020 Carlos M. Fonseca <cmfonsec@dei.uc.pt>
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
#include "shd.h"

extern gsl_rng *rng;    /* The single rng instance used by the whole code */

struct solverState {
    struct problem *p;
    struct solution *s;
    struct move *v;
};

struct solverState *newSolver(struct problem *p) {
    struct solverState *ss;
    if (getNumObjectives(p) != 1) {
        fprintf(stderr, "shd error: problem must have a single objective.");
        return NULL;
    }
    ss = (struct solverState *) malloc(sizeof (struct solverState));
    ss->p = p;
    ss->s = allocSolution(p);
    ss->v = allocMove(p);
    randomSolution(ss->s);
    return ss;
}

void freeSolver(struct solverState *ss) {
    freeSolution(ss->s);
    freeMove(ss->v);
    free(ss);
}

struct solverState *nextSolverState(struct solverState *ss) {
    double inc;
    randomMove(ss->v, ss->s);
    if (*getObjectiveIncrement(&inc, ss->v, ss->s) <= 0) {
        applyMove(ss->s, ss->v);
    }
    return ss;
}

struct solution *getSolverSolution(struct solverState *ss) {
    return ss->s;
}

void printSolverState(struct solverState *ss) {
    printSolution(ss->s);
}

