/* shd.c
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
#include <gsl/gsl_rng.h>
#include "shd.h"

extern gsl_rng *rng;    /* The single rng instance used by the whole code */

struct solverState {
    struct problem *p;
    struct solution *s0, *s1;
    struct move *v;
};

struct solverState *newSolver(struct problem *p) {
    struct solverState *ss;
    ss = (struct solverState *) malloc(sizeof (struct solverState));
    ss->p = p;
    ss->s0 = allocSolution(p);
    ss->s1 = allocSolution(p);
    ss->v = allocMove(p);
    randomSolution(ss->s0);
    /* Force initialisation of the whole state
    nextSolverState(ss); */
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
    copySolution(ss->s1, ss->s0);
    randomMove(ss->v, ss->s1);
    applyMove(ss->s1, ss->v);
    if (getObjectiveValue(ss->s1) <= getObjectiveValue(ss->s0)) {
        tmp = ss->s0;
        ss->s0 = ss->s1;
        ss->s1 = tmp;
    }
    return ss;
}

struct solution *getSolverSolution(struct solverState *ss) {
    return ss->s0;
}

void printSolverState(struct solverState *ss) {
    printSolution(ss->s0);
}

