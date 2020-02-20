/* nqueens.c
 *
 * (C) 2020 Carlos M. Fonseca <cmfonsec@dei.uc.pt> and
 *          Leonor Coelho <leonor.coelho@gmail.com> and
 *          Rúben Leal <rleal@student.dei.uc.pt> and
 *          Samuel Outeiro <souteiro@student.dei.uc.pt>
 * (C) 2018 Eva Tuba <etuba@ieee.org> and Carlos M. Fonseca <cmfonsec@dei.uc.pt>
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
#include <string.h>
#include "nqueens.h"

struct problem {
    float *table;
    int n;   /* number of queens */
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

struct segment {
    struct problem *prob;
    int *data;      /* normalised permutation computed as p2i[p1] */ 
    int *datai;     /* inverse permutation of data */ 
    int *cycle ;    /* number of the cycle each element belongs to */
    int *pos;       /* positions in which the solutions differ */
    int n;          /* number of queens / board size */
    int n_cycles;   /* number of cycles in data */
    int n_diff;     /* number of positions in which the solutions differ */
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
 * N-Queens instantiation
 * Status: TENTATIVE
 * Notes:
 *   Should just take an integer as an argument.
 *   Needs error checking
 */
struct problem *newProblem(int n) {
    struct problem *p = NULL;
    if (n > 0) {
        p = (struct problem *) malloc(sizeof (struct problem));
        p->n = n;
        p->table = (float *) malloc(n*(n+1)/2 * sizeof(float));
        init_table(p->table, n);
    } else
        fprintf(stderr, "nqueens: Invalid board size: %d\n", n);
#if 0
    int i, k;
    for (i = 0; i < n; i++) {
        for (k = 0; k <= i; k++)
            printf("%8.4f", p->table[ix(i,k)]);
        printf("\n");
    }
#endif
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
    s->n = n;
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
 * Allocate memory for a path
 * Status: TENTATIVE
 */
struct segment *allocSegment(struct problem *p) {
    struct segment *ps = malloc (sizeof (struct segment));
    ps->prob = p;
    ps->data = malloc(p->n * sizeof (int));
    ps->datai = malloc(p->n * sizeof (int));
    ps->cycle = malloc(p->n * sizeof (int));
    ps->pos = malloc(p->n * sizeof (int));
    ps->n = p->n;
    return ps;
}

/*
 * Free the memory used by a problem
 * Status: FINAL
 */
void freeProblem(struct problem *p) {
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
 * Free the memory used by a path
 * Status: IN_PROGRESS
 * Notes:
 *   update along with segment
 */
void freeSegment(struct segment *ps) {
    free(ps->data);
    free(ps->datai);
    free(ps->cycle);
    free(ps->pos);
    free(ps);
}

/* I/O  */
void printProblem(struct problem *p) {
    printf("%d-Queens problem\n",p->n);
}

void printSolution(struct solution *s) {
    int i;
    printf("%d-Queens solution\n  p:", s->n);
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
    printf("%d-Queens move: %d, %d\n", v->prob->n, v->data[0], v->data[1]);
}

void printSegment(struct segment *ps) {
    int i;
    printf("%d-Queens path state\n", ps->n);
    printf("  Path length: %d\n  Permutation:", ps->n - ps->n_cycles);
    for (i = 0; i < ps->n; i++)
        printf(" %d",ps->data[i]);
    printf("\n  Inverse permutation:");
    for (i = 0; i < ps->n; i++)
        printf(" %d", ps->datai[i]);
    printf("\n  Cycles:");
    for (i = 0; i < ps->n; i++)
        printf(" %d", ps->cycle[i]);
    printf("\n  Number of different positions: %d\n  Positions: ", ps->n_diff);
    for (i = 0; i < ps->n_diff; i++)
        printf("%d ",ps->pos[i]);
    printf("\n");
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
 * Solution evaluation
 * Status: INTERIM
 * Notes:
 *   Implements incremental evaluation for multiple moves
 */
double *getObjectiveVector(double *objv, struct solution *s) {
    int i, j, n = s->n, att, *mod, nmod = s->nmod;
    if (nmod > .28*n) { /* full evaluation threshold to be fine-tuned experimentally */
        att = 0;
        for (i = 0; i < n-1; i++)
            for (j = i+1; j < n; j++)
                att += (j - i == abs(s->data[i] - s->data[j]));
        s->objvalue = att;
        memcpy(s->olddata, s->data, n * sizeof (int));
        s->nmod = 0;
    } else if (nmod >= 1) { /* incremental evaluation */
        att = s->objvalue;
        mod = s->mod;
        for (i = 0; i < nmod-1; i++)
            for (j = i+1; j < nmod; j++)
                att += (abs(mod[j] - mod[i]) == abs(s->data[mod[i]] - s->data[mod[j]])) -
                    (abs(mod[j] - mod[i]) == abs(s->olddata[mod[i]] - s->olddata[mod[j]]));
        for (i = 0; i < nmod; i++)
            for (j = nmod; j < n; j++)
                att += (abs(mod[j] - mod[i]) == abs(s->data[mod[i]] - s->data[mod[j]])) -
                    (abs(mod[j] - mod[i]) == abs(s->olddata[mod[i]] - s->olddata[mod[j]]));
        s->objvalue = att;
        memcpy(s->olddata, s->data, n * sizeof (int));
        s->nmod = 0;
#if 0
        /* for testing */
        att = 0;
        for (i = 0; i < n-1; i++)
            for (j = i+1; j < n; j++)
                if (j - i == abs(s->data[i] - s->data[j]))
                    att++;
        printf("incremental = %d, full = %d, diff = %d\n", s->objvalue, att, s->objvalue-att);

#endif
    }
    *objv = (double)s->objvalue;
    return objv;
}

int getExcentricity(struct solution *s) {
    return s->n-1;
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
    int i, j, k, n = s->n, att;
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
    att = 0;
    /*
        It should be noted that:
        (abs(i - j) == abs(s->data[i] - s->data[j])) == (abs(i - j) == abs(s->data[j] - s->data[i])),
        that is, the two swapped queens interact the same way before and after swapping. Therefore,
        there is no need to compute their interaction.
    */
    for (k = 0; k < i; k++) {
        att += ((i - k) == abs(s->data[k] - s->data[j])) - ((i - k) == abs(s->data[k] - s->data[i]));
        att += ((j - k) == abs(s->data[k] - s->data[i])) - ((j - k) == abs(s->data[k] - s->data[j]));
    }
    for (k = i+1; k < j; k++) {
        att += ((k - i) == abs(s->data[k] - s->data[j])) - ((k - i) == abs(s->data[k] - s->data[i]));
        att += ((j - k) == abs(s->data[k] - s->data[i])) - ((j - k) == abs(s->data[k] - s->data[j]));
    }
    for (k = j+1; k < n; k++) {
        att += ((k - i) == abs(s->data[k] - s->data[j])) - ((k - i) == abs(s->data[k] - s->data[i]));
        att += ((k - j) == abs(s->data[k] - s->data[i])) - ((k - j) == abs(s->data[k] - s->data[j]));
    }
    *obji = (double)(v->incrvalue = att);
    return obji;
}

/*
 * Generate a random move in towards the second extreme of a segment
 * Status: CHECK, NEEDS_WORK
 * Notes:
 *  FIXME: should compute moves by using cycles
 */
struct move *randomMoveTowards(struct move *v, struct segment *ps) {
    int i, j, r, n = ps->n;
    int *pos = ps->pos;

    if (ps->n_cycles >= n) /* end of path has been reached, no move is possible */
        return NULL;

    /* Note: a swap move may correct two positions at once, but we only check
       one at a time. Therefore, we may need to skip null moves here, and try
       again. FIXME: input arg segment cannot be const for this reason! */
    do {
        /* Draw at random a position to correct */
        r = randint(--ps->n_diff);
        i = v->data[0] = pos[r];  /* choose an element that is not correct */ 
        /* Find where correct one is */
        j = ps->datai[i];
        pos[r] = pos[ps->n_diff];
    } while (i == j);
    v->data[1] = j;
    return v;
}

struct segment *applyMoveToSegment(struct segment *ps, const struct move *v) {
    int i, j;
    /* Update segment */
    i = v->data[0];
    j = v->data[1];
    /* FIXME: double check this, it can probably be written more simply */
    swap(ps->datai, ps->data[i], ps->data[j]);
    swap(ps->data, i, j);
    ps->n_cycles += 2 * (ps->cycle[i] == ps->cycle[j]) - 1;
    /* FIXME: Must update cycle[] here to reflect the move applied. */
    return ps;
}

/* Path generation */

/*
 * Set up a path from one solution to another
 * Status: CHECK
 * Notes:
 *   find and save cycles of permutation p2i[p1]
 */
struct segment *initSegment(struct segment *ps, const struct solution *s1, const struct solution *s2) {
    int i, j, n = ps->n;
    
    /* Compute reference permutation p2i[p1] using ps->pos as a buffer */
    for (i = 0; i < n; i++)
        ps->pos[s2->data[i]] = i;
    for (i = 0; i < n; i++)
        ps->data[i] = ps->pos[s1->data[i]];
    /* Inverse of reference permutation */
    for (i = 0; i < n; i++)
        ps->datai[ps->data[i]] = i;

    /* Compute the cycles in reference permutation p2i[p1] */
    ps->n_cycles = 0;
    for (i = 0; i < n; i++)
        ps->cycle[i] = 0;    /* all elements are unchecked */
    for (i = 0; i < n; i++) {
        if (ps->cycle[i] == 0) {
            ps->n_cycles++;
            j = i;
            do {
                ps->cycle[j] = ps->n_cycles;
                j = ps->data[j];
            } while (j != i);
        }
    }

    /* Find elements of data that differ from their position */
    for (i = 0, j = 0; i < n; i++)
        if (ps->data[i] != i)
            ps->pos[j++] = i;
    ps->n_diff = j;
    return ps;    
}

/* Path inspection */

/*
 * Current length of path
 * Status: FINAL
 */
int getLength(struct segment *ps) {
    return ps->n - ps->n_cycles;
}

