
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "symmtsp.h"
#include <gsl/gsl_rng.h>



struct problem{
	double *distances;	/* matrix holding distances between cities */
	int n; 				/* number of cities */
};

struct solution{
	struct problem *prob;
	int *data;
	int n;
	int evaluated; /* flag which indicates if solution has evaluated objValue */
	double objValue;
    int *movePool; /* array which holds all possible moves for current solution*/
    int movePoolLim; /* number of unused moves from the pool*/
};

struct move{
	struct problem *prob;
	int data[2];   		/* two positions of current permutation that will change place */
};


extern gsl_rng *rng;    /* The single rng instance used by the whole code */

static int randint(int n_max) {
    return gsl_rng_uniform_int(rng, n_max+1);
}

/**********************************/
/* ----- Utility functions ----- */
/**********************************/

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
    int el = data[i];
    data[i] = data[j];
    data[j] = el;
}

static int nh_size(int n) {
    return (n <= 3) ? 0 : n * (n - 3);
}



/*************************/
/* Problem instantiation */
/*************************/

struct problem *newProblem(const char *filename){
		int n, i, j;
		FILE *inFile;
		struct problem *p = NULL;
		inFile = fopen(filename, "r");

		if(inFile){
			fscanf(inFile, "%d", &n);
			if(n > 1){
				p = malloc(sizeof(struct problem));
				p -> n = n;
				p -> distances = malloc(sizeof(double) * n * n);
				for(i = 0; i < n; ++i){
					for(j = 0; j < n; ++j){
						fscanf(inFile, "%lf", &p -> distances[i * n + j]);

						if(j < i && p->distances[i * n + j] != p->distances[j * n + i]){
							fprintf(stderr, " Invalid traveling salesman problem formulation from %s.\n Distances matrix is not symmetric", filename);
							fclose(inFile);
							return NULL;
						}

					}
				}
			}else
				fprintf(stderr, " Invalid traveling salesman problem formulation from %s.\n Number of cities must be greater than 1", filename);

			fclose(inFile);
		}else{
			fprintf(stderr, "Cannot open file %s\n", filename);
			return p;
		}
		return p;
}


/**********************/
/* Problem inspection */
/**********************/

int getNumObjectives(const struct problem *p){
	return 1;
}


/*********************/
/* Memory management */
/*********************/

struct solution *allocSolution(struct problem *p){
	int i, j, n = p->n, k = 0;
	struct solution *s = malloc( sizeof(*s) );
	s -> data = malloc( sizeof(int*) * p->n );
	s -> prob = p;
	s -> n = n;
	s -> evaluated = 0;
	/* randomMoveWOR() support */
	s->movePoolLim = nh_size(n);
	s->movePool = malloc(sizeof(int) * s->movePoolLim);
	for(i = 0; i < n; ++i )
		for(j = i + 2; j < i + 2 + n - 3; ++j)
			s->movePool[k++] = i * n + j % n;
	return s;
}


struct move *allocMove(struct problem *p){
	struct move *m = malloc( sizeof(*m) );
	m -> prob = p;
	return m;
}

void freeProblem(struct problem *p){
	free(p -> distances);
	free(p);
}


void freeSolution(struct solution *s){
	free(s-> data);
	free(s->movePool);
	free(s);
}


void freeMove(struct move *v){
	free(v);
}


/*******/
/* I/O */
/*******/

void printProblem(struct problem *p){
	int i,j,n;
	n = p -> n;
	printf("--------------------------\n");
	printf("Printing problem:\n");
	printf("Matrix of distances:\nIndex     ");
	for(i = 0; i < n ; ++i)
		printf("%d        ", i);
	printf("\n");
	for(i = 0; i < n; ++i){
		printf("%d      ", i);
		for(j = 0; j < n; ++j)
			printf("%.3lf  ",*(p->distances + i * n + j));
		printf("\n");
	}
	printf("--------------------------\n");
}


void printSolution(struct solution *s){
	int i;
	printf("Solution ");
	for(i = 0; i < s->n ; ++i)
		printf("|%d|", *(s->data + i));
	if(s->evaluated)
		printf("  Objective value:%f", s->objValue);
	printf("\n");
}


void printMove(struct move *v){
	printf("Move  [%d,%d]\n", v->data[0], v->data[1]);
};



/***********************/
/* Solution generation */
/***********************/
struct solution *randomSolution(struct solution *s){
	randperm(s->data,s->n);
	s->evaluated = 0;
	s->movePoolLim = nh_size(s->n);
	return s;
}

//TODO develop other randomNeighbour strategy
struct solution *randomNeighbour(struct solution *s, int d){
	int i;
	struct move *v = allocMove(s->prob);
	for (i = 0; i < d; i++)
	   applyMove(s, randomMove(v, s));
	return s;
}


/***********************/
/* Solution inspection */
/***********************/
double *getObjectiveVector(double *objv, struct solution *s){
	int i;
	if(s->evaluated){
		*objv = s->objValue;
	}else{
		*objv = 0;
		for(i = 0; i < s->n - 1; ++i)
			*objv += s->prob->distances[s->data[i] * s->n + s->data[i+1]];
		*objv += s->prob->distances[s->data[s->n - 1] * s->n + s->data[0]];
	}
	return objv;
}


int equalSolutions(struct solution *s1, struct solution *s2){
	int i, startIndex;
	int f1, f2;//flags
	if(s1->prob != s2->prob)
		return 0;

	for(i = 0; i < s1->n; ++i){
		if( s1->data[0] == s2->data[i]){
			startIndex = i;
			break;
		}
	}
	f1 = f2 = 1;
	for(i = 1; i < s1->n; ++i){
		if( s1->data[i] != s2->data[ (startIndex + i)%s2->n ] )
			f1 = 0;
		if( s1->data[i] != s2->data[ (s2->n + startIndex - i)%s2->n] )
			f2 = 0;
		if(f1 == 0 && f2 == 0)
			return 0;
	}
	return 1;
}


int getNeighbourhoodSize(struct solution *s){
	return nh_size(s->n);
}

int isSolEvaluated(struct solution *s){
	return s->evaluated;
}


/*******************/
/* Move generation */
/*******************/
struct move *randomMove(struct move *v, const struct solution *s){
	if(s->n <= 3){
		v->data[0] = 0;
		v->data[1] = 1;
	}else{
		v->data[0] = randint(s->n-1);
		v->data[1] = (v->data[0] + 2 + randint(s->n-4)) % s->n;
	}
    return v;
}


struct move *randomMoveWOR(struct move *v, struct solution *s){
	int r, move;
	if(s->movePoolLim <= 0)
		return NULL;
	r = randint(--s->movePoolLim);
	move = s->movePool[r];  // get random unused move from the pool
	//swap move to the used moves part
	s->movePool[r] = s->movePool[s->movePoolLim];
	s->movePool[s->movePoolLim] = move;
	// formula for generating the move [i,j] was: i * n + j % n
	v->data[0] = move / s->n; // i
	v->data[1] = move % s->n; // j
	return v;
}


struct solution *resetRandomMoveWOR(struct solution *s){
	s->movePoolLim = nh_size(s->n);
	return s;
}




/**************************/
/* Operations on solutions*/
/**************************/
struct solution *copySolution(struct solution *dest, const struct solution *src){
	dest->prob = src->prob;
	dest->n = src->n;
	memcpy(dest->data, src->data, src->n * sizeof (int));
	memcpy(dest->movePool, src->movePool, nh_size(src->n) *sizeof(int));
	dest->movePoolLim = src->movePoolLim;
	dest->evaluated = src->evaluated;
	dest->objValue = src->objValue;
	return dest;
}


struct solution *applyMove(struct solution *s, const struct move *v){
	int i,nSteps;
	if(v->data[0] < v->data[1]){
		nSteps = (v->data[1] - v->data[0]) / 2;
		for(i = 0; i < nSteps; ++i)
			swap(s->data, v->data[0] + i, v->data[1] - 1 - i);
	}else{
		nSteps = (s->n - ( v->data[0] - v->data[1]) ) / 2;
		for(i = 0;i < nSteps; ++i)
			swap(s->data, (v->data[0] + i) % s->n, (s->n + v->data[1] - 1 - i) % s->n);
	}
	s->evaluated = 0;
	s->movePoolLim = nh_size(s->n);
	return s;
}



struct solution *setObjValue(double objv, struct solution *s){
	s->objValue = objv;
	s->evaluated = 1;
	return s;
}


/*******************/
/* Move inspection */
/*******************/
double *getObjectiveIncrement(double *obji, struct move *v, struct solution *s){
	int c1,c2, c1n, c2n;
	c1 = s->data[ v->data[0] ];
	c2 = s->data[ (s->n + v->data[1] - 1)%s->n ];
	c1n = s->data[ (s->n + v->data[0] - 1)%s->n ]; //city next to v1
	c2n = s->data[ v->data[1] ]; //city next to v2
	*obji = - s->prob->distances[c1 * s->n + c1n] - s->prob->distances[c2 * s->n + c2n]
		    + s->prob->distances[c2 * s->n + c1n] + s->prob->distances[c1 * s->n + c2n];
	return obji;
}


