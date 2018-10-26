#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "knapsack/knapsack.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

gsl_rng *rng;    /* The single rng instance used by the whole code */

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

void ranking(double *fitness, double **aux, double *cost, double sp, int n) {
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

void sus(struct solution **o, struct solution **p, const double *f, int nsel, int n) {
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
void single_step_mutation(struct solution **pop, int n, double Pm, struct move *tmpmove) {
    int i;
    for (i = 0; i < n; i++)
        if (gsl_rng_uniform(rng) < Pm) {
            randomMove(tmpmove, pop[i]);
            applyMove(pop[i], tmpmove);
        }
}

/* Standard mutation (in place) */
void standard_mutation(struct solution **pop, int n, double Pm, struct move *tmpmove, struct pathState *tmppath) {
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
void geometric_recombination(struct solution **pop, int n, double Px, struct move *tmpmove, struct pathState *tmppath) {
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

int main(int argc, char **argv) {
    struct problem *p;
    struct solution **parent, **offspring, **tmppop;
    struct move *v;
    struct pathState *ps;
    double *cost, *fitness, mincost, **ranking_aux;
    int i, mu, gen, maxgen;
    if (argc < 4) {
        printf("Usage: simple_solver <problem_instance_filename> <popsize> <maxgen>\n");
        return 0;
    }

    /* Set up random number generation */
    rng = gsl_rng_alloc(gsl_rng_mt19937); 
    gsl_rng_set(rng, time(0));

    /* Problem instantiation and solver set up */
    p = newProblem(argv[1]);
    if (p == NULL) {
        fprintf(stderr, "Invalid problem instance\n");
        gsl_rng_free(rng);
        return 1;
    }
    mu = atoi(argv[2]);
    maxgen = atoi(argv[3]);
    v = allocMove(p);
    ps = allocPathState(p);
    parent = malloc(mu * sizeof (struct solution *));
    offspring = malloc(mu * sizeof (struct solution *));
    for (i = 0; i < mu; i++) {
        parent[i] = allocSolution(p);
        offspring[i] = allocSolution(p);
    }
    cost = malloc(mu * sizeof (double));
    fitness = malloc(mu * sizeof (double));
    ranking_aux = malloc(mu * sizeof (double *));

    /* Population initialisation */
    for (i = 0; i < mu; i++) {
        randomSolution(parent[i]);
    }
    mincost = INFINITY;
    for (gen = 0; gen < maxgen; gen++) {

        /* Individual evaluation */
        for (i = 0; i < mu; i++) {
            cost[i] = getObjectiveValue(parent[i]);
            if (cost[i] < mincost) {
                mincost = cost[i];
                printf("gen = %d, obj = %.0f\n", gen, mincost);
            }
        }

        /* Fitness assignment */
        ranking(fitness, ranking_aux, cost, 2, mu);

        /* Parental selection (sampling) */
        sus(offspring, parent, fitness, mu, mu);

        /* Geometric recombination (in place) */
        geometric_recombination(offspring, mu, 0.6, v, ps);

#if 1
        /* Single-step mutation (in place) */
        single_step_mutation(offspring, mu, (1 - 1./2.)*.8, v);
#else
        /* Independent mutation (in place) */
        standard_mutation(offspring, mu, (1 - 1./2.)*.8, v, ps);
#endif

        /* Unconditional generational replacement */
        tmppop = parent;
        parent = offspring;
        offspring = tmppop;
    }

    /* Free allocated memory */
    free(ranking_aux);
    free(fitness);
    free(cost);
    for (i = 0; i < mu; i++) {
        freeSolution(offspring[i]);
        freeSolution(parent[i]);
    }
    free(offspring);
    free(parent);
    freePathState(ps);
    freeMove(v);
    freeProblem(p);
    gsl_rng_free(rng);
    return 0;
}
