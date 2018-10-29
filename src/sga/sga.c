/* sga.c
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
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "sga.h"

extern gsl_rng *rng;    /* The single rng instance used by the whole code */

struct solverState {
    struct problem *p;
    struct solution **parent, **offspring, *best;
    struct move *v;
    struct pathState *ps;
    double *cost, *fitness, mincost, **ranking_aux;
    int popsize, i;
    double Px, Pm, SP;
};

/* Auxiliary functions */

static int cmp_doublep(const void *a, const void *b) {
    double x = **(double **) a;
    double y = **(double **) b;
    if (x < y)
        return -1;
    else if (x > y)
        return 1;
    else
        return 0;
}

static void ranking(double *fitness, double **aux, double *cost, double sp, int n) {
    int i, j, k;
    double num, den;
    for (i = 0; i < n; i++)
        aux[i] = cost + i;
    qsort(aux, n, sizeof(double *), cmp_doublep);
    for (i = 0; i < n; i = j) {
        num = sp - (sp-1.) * 2.0 * i / (n-1.);
        den = 1.;
        for (j = i+1; j < n && *aux[i] == *aux[j]; j++) {
            num += sp - (sp-1.) * 2.0 * j / (n-1.);
            den += 1.;
        }
        for (k = i; k < j; k++)
            fitness[aux[k]-cost] = num / den;
    }
}

static void sus(struct solution **o, struct solution **p, const double *f, int nsel, int n) {
    struct solution *tmp;
    double x, cumfit, totfit = 0.0;
    int i, j, last_i;

    for (i = j = 0; i < n; i++)
        totfit += f[i];
    x = gsl_rng_uniform(rng);
    cumfit = f[0];
    last_i = -1;
    i = j = 0;
    while (j < nsel) {
        if (totfit * x < cumfit * nsel) {
            if (i != last_i) { /* avoid copying the first time */
                tmp = o[j];
                o[j] = p[i];
                p[i] = tmp;
                last_i = i;
            } else
                copySolution(o[j], o[j-1]);
            j++;
            x += 1.;
        } else if (i < n-1) {
            cumfit += f[++i];
        } else {
            /* should never *ever* happen */
            fprintf(stderr, "SUS warning: end of parent array reached.");
            break;
        }
    }
    /* Shuffle offspring */
    for (i = n-1; i > 0; i--) {
        j = gsl_rng_uniform_int(rng, i+1);
        tmp = o[i];
        o[i] = o[j];
        o[j] = tmp;
    }
}

/* Single step mutation (in place) */
static void single_step_mutation(struct solution **pop, int n, double Pm, struct move *tmpmove) {
    int i;
    for (i = 0; i < n; i++)
        if (gsl_rng_uniform(rng) < Pm) {
            randomMove(tmpmove, pop[i]);
            applyMove(pop[i], tmpmove);
        }
}

/* Standard mutation (in place) */
static void standard_mutation(struct solution **pop, int n, double Pm, struct move *tmpmove, struct pathState *tmppath) {
    double pm;
    int i, j, k;
    for (i = 0; i < n; i++) {
        initPathAwayFrom(tmppath, pop[i]);
        k = getPathLength(tmppath);
        pm = 1.-pow(1.-Pm, 1./k);
        k = gsl_ran_binomial(rng, pm, k);
        for (j = 0; j < k; j++) {
            nextRandomMove(tmpmove, tmppath);
            applyMove(pop[i], tmpmove);
        }
    }
}

/* Geometric recombination (in place) */
static void geometric_recombination(struct solution **pop, int n, double Px, struct move *tmpmove, struct pathState *tmppath) {
    struct solution *tmpsol;
    int i, k;
    tmpsol = pop[n-1];
    for (i = 0; i < n; i++) {
        if (gsl_rng_uniform(rng) < Px) {
            initPathTo(tmppath, tmpsol, pop[i]);
            k = getPathLength(tmppath);
            if (k > 1) {
                k -= gsl_rng_uniform_int(rng, k/2)+1;
                /* This depends on the current path length being updated by
                 * nextRandomMove() and assumes that calling getPathLength() is
                 * cheap. */
                while (getPathLength(tmppath) > k) {
                    nextRandomMove(tmpmove, tmppath);
                    applyMove(tmpsol, tmpmove);
                }
            }
        }
        tmpsol = pop[i];
    }
}

/* Solver API functions */

struct solverState *newSolver(struct problem *p, int popsize) {
    struct solverState *ss;
    int i;
    /* memory allocation */
    ss = malloc(sizeof (struct solverState));
    ss->p = p;
    ss->popsize = popsize;
    ss->parent = malloc(popsize * sizeof (struct solution *));
    for (i = 0; i < popsize; i++)
        ss->parent[i] = allocSolution(p);
    ss->best = allocSolution(p);
    ss->cost = malloc(popsize * sizeof (double));
    /* temporary workspace */
    ss->offspring = malloc(popsize * sizeof (struct solution *));
    for (i = 0; i < popsize; i++)
        ss->offspring[i] = allocSolution(p);
    ss->v = allocMove(p);
    ss->ps = allocPathState(p);
    ss->ranking_aux = malloc(popsize * sizeof (double *));
    ss->fitness = malloc(popsize * sizeof (double));
    /* state initialisation */
    for (i = 0; i < popsize; i++) {
        randomSolution(ss->parent[i]);
    }
    ss->mincost = INFINITY;
    for (i = 0; i < popsize; i++) {
        ss->cost[i] = getObjectiveValue(ss->parent[i]);
        if (ss->cost[i] < ss->mincost) {
            ss->mincost = ss->cost[i];
            copySolution(ss->best, ss->parent[i]);
        }
    }
    /* default parameter values */
    ss->SP = 2.;
    ss->Px = .6;
    ss->Pm = (1.-1./ss->SP)*.8;
    return ss;
}

void freeSolver(struct solverState *ss) {
    int i;
    free(ss->ranking_aux);
    free(ss->fitness);
    freePathState(ss->ps);
    freeMove(ss->v);
    for (i = 0; i < ss->popsize; i++) {
        freeSolution(ss->offspring[i]);
        freeSolution(ss->parent[i]);
    }
    freeSolution(ss->best);
    free(ss->cost);
    free(ss->offspring);
    free(ss->parent);
    free(ss);
}

struct solverState *nextSolverState(struct solverState *ss) {
    struct solution **tmppop;
    int popsize = ss->popsize, i;

    /* Fitness assignment */
    ranking(ss->fitness, ss->ranking_aux, ss->cost, ss->SP, ss->popsize);

    /* Parental selection (sampling) */
    sus(ss->offspring, ss->parent, ss->fitness, popsize, popsize);

    /* Geometric recombination (in place) */
    geometric_recombination(ss->offspring, popsize, ss->Px, ss->v, ss->ps);

#if 1
    /* Single-step mutation (in place) */
    single_step_mutation(ss->offspring, popsize, ss->Pm, ss->v);
#else
    /* Independent mutation (in place) */
    standard_mutation(ss->offspring, popsize, ss->Pm, ss->v, ss->ps);
#endif

    /* Unconditional generational replacement */
    tmppop = ss->parent;
    ss->parent = ss->offspring;
    ss->offspring = tmppop;
    for (i = 0; i < popsize; i++) {
        ss->cost[i] = getObjectiveValue(ss->parent[i]);
        if (ss->cost[i] < ss->mincost) {
            ss->mincost = ss->cost[i];
            copySolution(ss->best, ss->parent[i]);
        }
    }
    return ss;
}

struct solution *getSolverSolution(struct solverState *ss) {
    return ss->best;
}

void printSolverState(struct solverState *ss) {
    fprintf(stderr, "sga: printSolverState() not implemented yet.\n");
}

/* Solver specific functions */

double getSelectivePressure(struct solverState *ss) {
    return ss->SP;
}

double getRecombinationRate(struct solverState *ss) {
    return ss->Px;
}

double getMutationRate(struct solverState *ss) {
    return ss->Pm;
}

void setSelectivePressure(struct solverState *ss, double SP) {
    if (SP >= 1. && SP <= 2.)
        ss->SP = SP;
    else
        fprintf(stderr, "sga: Invalid selective pressure. Unchanged.\n");
}

void setRecombinationRate(struct solverState *ss, double Px) {
    if (Px >= 0. && Px <= 1.)
        ss->Px = Px;
    else
        fprintf(stderr, "sga: Invalid recombination rate. Unchanged.\n");
}

void setMutationRate(struct solverState *ss, double Pm) {
    if (Pm >= 0. && Pm <= 1.)
        ss->Pm = Pm;
    else
        fprintf(stderr, "sga: Invalid mutation rate. Unchanged.\n");
}

