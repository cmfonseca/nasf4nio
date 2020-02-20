/* qap.c
 *
 * (C) 2020 Carlos M. Fonseca <cmfonsec@dei.uc.pt> and
 *          Leonor Coelho <leonor.coelho@gmail.com> and
 *          Rúben Leal <rleal@student.dei.uc.pt> and
 *          Samuel Outeiro <souteiro@student.dei.uc.pt>
 *
 * Modelled after nqueens.c
 * (C) 2018 Eva Tuba <etuba@ieee.org> and Carlos M. Fonseca <cmfonsec@dei.uc.pt>
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
#include <stdio.h>
#include <string.h>
#include "qap.h"

struct problem {
    int *flow, *dist;
    float *table;
    int n;   /* problem size */
};

struct solution {
    struct problem *prob;
    int *data;
    int *olddata;
    int *mod;       /* positions modified by moves */
    int *modi;      /* inverse permutation of mod */
    int nmod;       /* no. of positions modified since last evaluation */
    int n;
    int objvalue;
    int *rndSample; /* array used for sampling without replacement */
    int sampleLim; 
};

struct move {
    struct problem *prob;
    int data[2];
    int incrvalue;
};

/**********************************/
/* ----- Utility functions ----- */
/**********************************/

#if 0
/*
 * Random integer x such that 0 <= x <= n_max
 * Status: FINAL
 */
static int randint(int n_max) {
    int x;
    if (n_max < 1)
        return 0;
    div_t y = div(-RAND_MAX-1, -n_max-1);
    do
        x = random();
    while (x > RAND_MAX + y.rem);
    return x / y.quot;
}
#else
#include <gsl/gsl_rng.h>
extern gsl_rng *rng;    /* The single rng instance used by the whole code */

static int randint(int n_max) {
    return gsl_rng_uniform_int(rng, n_max+1);
}
#endif

/*
 * Random permutation of size n
 * Status: FINAL
 * Notes:
 *   Based on https://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle#The_%22inside-out%22_algorithm
 */
static int *randperm(int *p, int n) {
    /* Inside-out algorithm */
    int i, j;
    if (n > 0) {
        p[0] = 0;
        for (i = 1; i < n; i++) {
            j = randint(i);
            p[i] = p[j];
            p[j] = i; /* = source[i] */
        }
    }
    return p;
}

static void swap(int *data, int i, int j) {
    if (i == j)
        return;
    int el = data[i];
    data[i] = data[j];
    data[j] = el;
}

static int nh_size(int n) {
    return n*(n-1)/2;
}

/*
 * S(n,k) are the Stirling numbers
 * B(n,k) = S(n+1,k)/S(n+1,k+1) = B(n-1,k)[n + B(n-1,k-1)]/[n + B(n-1,k)]
 *
 * P(n,k) = S(n,k)/S(n+1,k+1) = B(n-1,k)/[n + B(n-1,k)]
 *
 * Create a unidimensional vector. First fill it with B numbers. Then,
 * start from the end to the beginning to calculate the probabilities.
*/

static int ix(int n, int k) {
    return n*(n+1)/2+k;
}

static void init_table(float *t, int n) {
    int i, k;
    t[0] = 0.;
    t[1] = 0.;
    t[2] = 1.;
    for(i = 2; i < n-1; i++) {
        t[ix(i,0)] = 0.;
        for(k = 1; k < i; k++)
            t[ix(i,k)] = t[ix(i-1,k)] *
                    (i + t[ix(i-1,k-1)]) / (i + t[ix(i-1,k)]);
        t[ix(i,i)] = i*(i+1)/2;
    }
    t[ix(n-1,0)] = 0.;
    for(i = n-1; i > 1; i--) {
        for(k = 1; k < i; k++)
            t[ix(i,k)] = t[ix(i-1,k)] / (i + t[ix(i-1,k)]);
        t[ix(i,i)] = 1.;
    }
    t[0] = 1.;
}

/**********************************************/
/* ----- Problem-specific instantiation ----- */
/**********************************************/

/*
 * QAP instantiation
 * Status: TENTATIVE
 * Notes:
 *   Needs error checking
 */
struct problem *newProblem(const char *filename) {
    struct problem *p = NULL;
    FILE *infile;
    int j, i, n = -1;
    infile = fopen(filename, "r");
    if (infile) {
        fscanf(infile, "%d", &n);
        if (n > 0) {
            p = malloc(sizeof (struct problem));
            p->flow = (int *) malloc(n * n * sizeof (int));
            p->dist = (int *) malloc(n * n * sizeof (int));
            for (i = 0; i < n; i++)
                for (j = 0; j < n; j++)
                    fscanf(infile, "%d", p->flow + n * i + j);
            for (i = 0; i < n; i++)
                for (j = 0; j < n; j++)
                    fscanf(infile, "%d", p->dist + n * i + j);
            p->n = n;
            p->table = (float *) malloc(n*(n+1)/2 * sizeof(float));
            init_table(p->table, n);
        } else
            fprintf(stderr, "Invalid QAP instance %s\n", filename);
        fclose(infile);
    } else
        fprintf(stderr, "Cannot open file %s\n", filename);
    return p;
}

int getNumObjectives(const struct problem *p) {
    return 1;
}

/*****************************/
/* ----- API functions ----- */
/*****************************/

/* Memory management */

/*
 * Allocate memory for a solution
 * Status: CHECK
 */
struct solution *allocSolution(struct problem *p) {
    int i, j, k = 0, n = p->n;
    struct solution *s = malloc(sizeof (struct solution));
    s->prob = p;
    s->data = malloc(n * sizeof (int));
    s->olddata = malloc(n * sizeof (int));
    s->mod = malloc(n * sizeof (int));
    s->modi = malloc(n * sizeof (int));
    s->n = p->n;
    /* randomMoveWOR() support */
    s->sampleLim = nh_size(n);
    s->rndSample = malloc(s->sampleLim * sizeof (int));
    for(j = 1; j < n; j++)
        for(i = 0; i < j; i++)
            s->rndSample[k++] = j*n+i;

    return s;
}

/*
 * Allocate memory for a move
 * Status: FINAL
 */
struct move *allocMove(struct problem *p) {
    struct move *v = malloc(sizeof (struct move));
    v->prob = p;
    return v;
}

/*
 * Free the memory used by a problem
 * Status: FINAL
 */
void freeProblem(struct problem *p) {
    free(p->dist);
    free(p->flow);
    free(p->table);
    free(p);
}

/*
 * Free the memory used by a solution
 * Status: FINAL
 */
void freeSolution(struct solution *s) {
    free(s->data);
    free(s->olddata);
    free(s->mod);
    free(s->modi);
    free(s->rndSample);
    free(s);
}

/*
 * Free the memory used by a move
 * Status: FINAL
 */
void freeMove(struct move *v) {
    free(v);
}

/*

/* I/O  */
void printProblem(struct problem *p) {
    int i, j;
    int n = p->n;
    printf("# Quadratic Assignment Problem (size %d)\n", n);
    printf("%d\n\n", n);
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++)
            printf(" %d",p->flow[i * n + j]);
        printf("\n");
    }
    printf("\n");
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++)
            printf(" %d",p->dist[i * n + j]);
        printf("\n");
    }
}

void printSolution(struct solution *s) {
    int i;
    printf("# QAP solution of size %d\n  p:", s->n);
    for(i = 0; i < s->n; i++)
        printf(" %d", s->data[i]);
#if 0
    printf("\n  mod :");
    for(i = 0; i < s->n; i++)
        printf(" %d", s->mod[i]);
    printf("\n  modi:");
    for(i = 0; i < s->n; i++)
        printf(" %d", s->modi[i]);
#endif
    printf("\n");
}

void printMove(struct move *v) {
    printf("QAP (%d) move: %d, %d\n", v->prob->n, v->data[0], v->data[1]);
}

/* Solution generation */

/*
 * Generate solutions uniformly at random
 * Status: CHECK
 */
struct solution *randomSolution(struct solution *s) {
    /* solution s must have been allocated with allocSolution() */
    int i;
    randperm(s->data, s->n);
    /* Solution not evaluated yet */
    for (i = 0; i < s->n; i++)
        s->mod[i] = s->modi[i] = i;
    s->nmod = s->n;
    s->sampleLim = nh_size(s->n);
    return s;
}

struct solution *randomNeighbour(struct solution *s, int d){
    float *table = s->prob->table;
    int j, k = s->n-d;
    for(int i = s->n-1; i >= k /* && k > 0 */; i--){
        if(gsl_rng_uniform(rng) >=  table[ix(i,k-1)]) {
            j = randint(i-1);
            swap(s->data, i, j);
            if (s->modi[i] >= s->nmod) {
                swap(s->mod, s->modi[i], s->nmod);
                swap(s->modi, i, s->mod[s->modi[i]]);
                s->nmod++;
            }
            if (s->modi[j] >= s->nmod) {
                swap(s->mod, s->modi[j], s->nmod);
                swap(s->modi, j, s->mod[s->modi[j]]);
                s->nmod++;
            }
        }
        else
            k -= 1;
    }
    s->sampleLim = nh_size(s->n);
    return s;
}

/* Solution inspection */

/*
 * Produce the objective value of a solution
 * Status: INTERIM
 * Notes:
 *   Implements incremental evaluation for multiple moves
 */
double *getObjectiveVector(double *objv, struct solution *s) {
    int i, j, n = s->n, obj, *mod, nmod = s->nmod;
    int *flow = s->prob->flow, *dist = s->prob->dist;
    if (nmod > (1 - 0.707) * n) { /* full evaluation threshold to be fine-tuned experimentally */
        obj = 0;
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                obj += flow[n * i + j] * dist[n * s->data[i] + s->data[j]];
        s->objvalue = obj;
        memcpy(s->olddata, s->data, n * sizeof (int));
        s->nmod = 0;
    } else if (nmod >= 1) { /* incremental evaluation */
        obj = s->objvalue;
        mod = s->mod;
        for (i = 0; i < nmod; i++)
            for (j = 0; j < n; j++)
                obj += flow[n * mod[i] + j] * (dist[n * s->data[mod[i]] + s->data[j]] -
                       dist[n * s->olddata[mod[i]] + s->olddata[j]]);
        for (i = nmod; i < n; i++)
            for (j = 0; j < nmod; j++)
                obj += flow[n * mod[i] + mod[j]] * (dist[n * s->data[mod[i]] + s->data[mod[j]]] -
                       dist[n * s->olddata[mod[i]] + s->olddata[mod[j]]]);
        s->objvalue = obj;
        memcpy(s->olddata, s->data, n * sizeof (int));
        s->nmod = 0;
#if 0
        /* for testing */
        obj = 0;
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                obj += flow[n * i + j] * dist[n * s->data[i] + s->data[j]];
        printf("incremental = %d, full = %d, diff = %d\n", s->objvalue, obj, s->objvalue-obj);
#endif
    }
    *objv = (double)s->objvalue;
    return objv;
}

int getNeighbourhoodSize(struct solution *s){
    return nh_size(s->n);
}

/* Operations on solutions*/
struct solution *copySolution(struct solution *dest, const struct solution *src) {
    dest->prob = src->prob;
    dest->n = src->n;
    memcpy(dest->data, src->data, src->n * sizeof (int));
    memcpy(dest->olddata, src->olddata, src->n * sizeof (int));
    memcpy(dest->mod, src->mod, src->n * sizeof (int));
    memcpy(dest->modi, src->modi, src->n * sizeof (int));
    memcpy(dest->rndSample, src->rndSample, nh_size(src->n) * sizeof (int));
    dest->sampleLim = src->sampleLim;
    dest->nmod = src->nmod;
    dest->objvalue = src->objvalue;
    return dest;
}

/*
 * Apply a move to a solution
 * Status: FINAL
 */
struct solution *applyMove(struct solution *s, const struct move *v) {
    int i, j;
    i = v->data[0];
    j = v->data[1];
    if (i == j)     /* do nothing */
        return s;
    swap(s->data, i, j);
    if (s->modi[i] >= s->nmod) {
        swap(s->mod, s->modi[i], s->nmod);
        swap(s->modi, i, s->mod[s->modi[i]]);
        s->nmod++;
    }
    if (s->modi[j] >= s->nmod) {
        swap(s->mod, s->modi[j], s->nmod);
        swap(s->modi, j, s->mod[s->modi[j]]);
        s->nmod++;
    }
    s->sampleLim = nh_size(s->n);
    return s;
}

/* Move generation */

/*
 * Generate moves uniformly at random
 * Status: TENTATIVE, NEEDS_TESTING
 * Notes:
 *   Move (i,j) such that i != j. Order is irrelevant and not enforced.
 */
struct move *randomMove(struct move *v, const struct solution *s) {
    /* move v must have been allocated with allocMove() */
    int n;
    n = s->n;
    v->data[0] = randint(n-1);
    v->data[1] = (v->data[0] + 1 + randint(n-2)) % n;
    return v;
}

struct solution *resetRandomMoveWOR(struct solution *s){
    s->sampleLim = nh_size(s->n);
    return s;
}

struct move *randomMoveWOR(struct move *v, struct solution *s) {
    /* move v must have been allocated with allocMove() */
    int r, x, n = s->n;
    //printf("Sample Lim: %d\r", s->sampleLim);
    if (s->sampleLim <= 0) return NULL;
    r = randint(--s->sampleLim);
    x = s->rndSample[r];
    s->rndSample[r] = s->rndSample[s->sampleLim];
    s->rndSample[s->sampleLim] = x;
    v->data[0] = x % n;
    v->data[1] = x / n;
    return v;
}

double *getObjectiveIncrement(double *obji, struct move *v, struct solution *s) {
    int i, j, k, n = s->n, obj;
    int *flow = s->prob->flow, *dist = s->prob->dist;
    if (v->data[0] == v->data[1]) { /* do nothing */
        *obji = 0.0;
        return obji;
    } else if (v->data[0] < v->data[1]) {   
        i = v->data[0];
        j = v->data[1];
    } else {
        i = v->data[1];
        j = v->data[0];
    }
    obj = 0;
    obj += flow[n * i + j] * (dist[n * s->data[j] + s->data[i]] - dist[n * s->data[i] + s->data[j]]);
    obj += flow[n * j + i] * (dist[n * s->data[i] + s->data[j]] - dist[n * s->data[j] + s->data[i]]);
    /* The following two lines should have no effect, unless flow[n*k1+k1] and/or dist[n*k2+k2] are not zero for some k1,k2 */
    obj += flow[n * i + i] * (dist[n * s->data[j] + s->data[j]] - dist[n * s->data[i] + s->data[i]]);
    obj += flow[n * j + j] * (dist[n * s->data[i] + s->data[i]] - dist[n * s->data[j] + s->data[j]]);
    for (k = 0; k < i; k++) {
        obj += flow[n * i + k] * (dist[n * s->data[j] + s->data[k]] - dist[n * s->data[i] + s->data[k]]);
        obj += flow[n * j + k] * (dist[n * s->data[i] + s->data[k]] - dist[n * s->data[j] + s->data[k]]);
        obj += flow[n * k + i] * (dist[n * s->data[k] + s->data[j]] - dist[n * s->data[k] + s->data[i]]);
        obj += flow[n * k + j] * (dist[n * s->data[k] + s->data[i]] - dist[n * s->data[k] + s->data[j]]);
    }
    for (k = i+1; k < j; k++) {
        obj += flow[n * i + k] * (dist[n * s->data[j] + s->data[k]] - dist[n * s->data[i] + s->data[k]]);
        obj += flow[n * j + k] * (dist[n * s->data[i] + s->data[k]] - dist[n * s->data[j] + s->data[k]]);
        obj += flow[n * k + i] * (dist[n * s->data[k] + s->data[j]] - dist[n * s->data[k] + s->data[i]]);
        obj += flow[n * k + j] * (dist[n * s->data[k] + s->data[i]] - dist[n * s->data[k] + s->data[j]]);
    }
    for (k = j+1; k < n; k++) {
        obj += flow[n * i + k] * (dist[n * s->data[j] + s->data[k]] - dist[n * s->data[i] + s->data[k]]);
        obj += flow[n * j + k] * (dist[n * s->data[i] + s->data[k]] - dist[n * s->data[j] + s->data[k]]);
        obj += flow[n * k + i] * (dist[n * s->data[k] + s->data[j]] - dist[n * s->data[k] + s->data[i]]);
        obj += flow[n * k + j] * (dist[n * s->data[k] + s->data[i]] - dist[n * s->data[k] + s->data[j]]);
    }
    /* for testing */
    // printf("%d\n", obj);
    *obji = (double)(v->incrvalue = obj);
    return obji;
}
