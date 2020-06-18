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
	int evaluated;		/* flag which indicates if solution has evaluated objValue */
	double objValue;
    int *movePool; 		/* array which holds all possible moves for current solution*/
    int movePoolLim; 	/* number of unused moves from the pool*/
};

struct move{
	struct problem *prob;
	int data[2];   		/* two positions of current permutation that will change place */
};


struct segment {
    struct problem *prob;
    int n;
    int *data;			/* difference in elements position of one solution to second solution */
    int *datai;			/* inverse of difference */
	int *breakPoints;	/* array which holds position of active break points */
    int *breakPointsi;	/* inverse of breakPoints */
	int *bpAtPosition;  /* is there break point at given position */
    int numBp;          /* number of active break points */
};


extern gsl_rng *rng;    /* The single rng instance used by the whole code */


/**********************************/
/* ----- Utility functions ----- */
/**********************************/

static int randint(int n_max) {
    return gsl_rng_uniform_int(rng, n_max+1);
}


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

static int *reversal(int *data, struct move *v, int n){
	int nSteps, i;
	if(v->data[0] < v->data[1]){
		nSteps = (v->data[1] - v->data[0]) / 2;
		for(i = 0; i < nSteps; ++i)
			swap(data, v->data[0] + i, v->data[1] - 1 - i);
	}else{
		nSteps = (n - ( v->data[0] - v->data[1]) ) / 2;
		for(i = 0;i < nSteps; ++i)
			swap(data, (v->data[0] + i) % n, (n + v->data[1] - 1 - i) % n);
	}
	return data;
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


struct segment *allocSegment(struct problem *p){
	struct segment *seg = malloc(sizeof(*seg));
	seg->prob = p;
	seg->data = malloc(sizeof(int) * p->n);
	seg->datai = malloc(sizeof(int) * p->n);
	seg->breakPoints = malloc(sizeof(int) * p->n);
	seg->breakPointsi = malloc(sizeof(int) * p->n);
	seg->bpAtPosition = malloc(sizeof(int) * p->n);
	seg->numBp = 0;
	seg->n = p->n;
	return seg;
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


void freeSegment(struct segment *seg){
	free(seg->data);
	free(seg->datai);
	free(seg->breakPoints);
	free(seg->breakPointsi);
	free(seg->bpAtPosition);
	free(seg);
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


void printSegment(struct segment *seg){
	int i;
	printf("-----------------------------\n");
	printf("Segment\n");
	printf("Data: 					");
	for(i = 0; i < seg->n; ++i)
		printf("%d ", seg->data[i]);
	printf("\n");
	printf("Number of active breakpoints:		%d\n", seg->numBp);
	printf("Positions of active breakpoints:	");
	for(i = 0; i < seg->numBp; ++i)
		printf("%d ", seg->breakPoints[i]);
	printf("\n-----------------------------\n");
}

/***********************/
/* Solution generation */
/***********************/
struct solution *randomSolution(struct solution *s){
	randperm(s->data,s->n);
	s->evaluated = 0;
	s->movePoolLim = nh_size(s->n);
	return s;
}


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
	reversal(s->data, v, s->n);
	s->evaluated = 0;
	s->movePoolLim = nh_size(s->n);
	return s;
}



struct solution *setObjValue(double objv, struct solution *s){
	s->objValue = objv;
	s->evaluated = 1;
	return s;
}


/**********************/
/* Segment operations */
/**********************/

struct segment *initSegment(struct segment *seg, const struct solution *s1, const struct solution *s2){
	int i,n = s1->n, diff;
	int *s2_datainverse;

	s2_datainverse = malloc(sizeof(int) * n);
	for(i = 0; i < n; ++i)
		s2_datainverse[s2->data[i]] = i;

	for(i = 0; i < n; ++i){
		seg->data[i] = s2_datainverse[s1->data[i]];
		seg->datai[seg->data[i]] = i;
	}

	seg->numBp = 0;
	// Check if elements at first and last position are adjacent
	diff =  abs( seg->data[0] - seg->data[n - 1] );
	if(diff != 1 && diff != (n - 1) ){
		seg->breakPoints[seg->numBp++] = 0;
		seg->bpAtPosition[0] = 1;
		seg->breakPointsi[0] = 0;
	}else{
		seg->bpAtPosition[0] = 0;
		seg->breakPointsi[0] = -1;
	}

	for(i = 1; i < n; ++i){
		diff = abs( seg->data[i] - seg->data[i - 1] );
		if( diff != 1 && diff != (n - 1) ){
			seg->breakPoints[seg->numBp] = i;
			seg->bpAtPosition[i] = 1;
			seg->breakPointsi[i] = seg->numBp;
			seg->numBp ++;
		}
		else{
			seg->bpAtPosition[i] = 0;
			seg->breakPointsi[i] = -1;
		}
	}
	free(s2_datainverse);
	return seg;
}


/*
 * Function checks if there is adjacent element to one at bp, smaller
 * or bigger, whose position is also on the break point.
 * If there is such an element, position of that element is returned.
 * If there is no such element, -1 is returned.
 */
int findSuitableBreakPoint(struct segment *seg, int bp){
	int elemAtBp, biggerAdjPos, smallerAdjPos , n = seg->n;

	elemAtBp = seg->data[bp];
	biggerAdjPos = seg->datai[(elemAtBp + 1) % n];
	smallerAdjPos = seg->datai[(elemAtBp - 1 + n) % n];

	if(seg->bpAtPosition[biggerAdjPos] )
		return biggerAdjPos;
	else if(seg->bpAtPosition[smallerAdjPos] )
		return smallerAdjPos;
	else
		return -1;
}

/*
 * Function randomly chooses breakpoints from the set of active breakpoints
 * and tries to find second suitable breakpoint which is also from the set of
 * active breakpoints.
 * If no suitable breakpoint is found, two randomly chosen breakpoints are set
 * for a move.
 *
 * Two breakpoints are suitable for a move if elements at those breakpoints are
 * adjacent. If they are, applying reversal will lead to reduction of at least
 * one breakpoint.
 */
struct move *randomMoveTowards(struct move *v, struct segment *seg){
	int r,bp1,bp2, n = seg->n, numBpCopy = seg->numBp;
	int *breakPointsCopy = malloc(sizeof(int) * seg->numBp);
	memcpy(breakPointsCopy ,seg->breakPoints, seg->numBp * sizeof(int));

	if(seg->numBp == 0)
		return NULL;


	while(numBpCopy){
		r = randint(numBpCopy - 1);
		bp1 = breakPointsCopy[r];
		bp2 = findSuitableBreakPoint(seg, bp1);
		if(bp2 != -1){
			v->data[0] = bp1;
			v->data[1] = bp2;
			free(breakPointsCopy);
			return v;
		}
		swap(breakPointsCopy, r, --numBpCopy);
	}
	// If non was found pick 2 randomly
	r = randint(seg->numBp - 1);
	v->data[0] = seg->breakPoints[r];
	v->data[1] = seg->breakPoints[ (r + 1 + randint(seg->numBp - 2)) % seg->numBp ];
	free(breakPointsCopy);
	return v;
}


struct segment *applyMoveToSegment(struct segment *seg, const struct move *v){
	int i, j, diff, n = seg->n;
	// Applying the move to data
	reversal(seg->data, v, seg->n);

	// Apply change that occurred to datai, bpAtPosition, breakPoints and breakPointsi
	for(i = (v->data[0] - 1 + n) % n; i != v->data[1]; i = j ){
		j = (i + 1) % n;
		seg->datai[seg->data[j]] = j;

		diff = abs(seg->data[i] - seg->data[j]);
		if( diff == 1 || diff == (n - 1) ){ // If elements are adjacent
			if( seg->bpAtPosition[j] ){  	// And there is leftover breakpoint in between them, remove breakpoint
				seg->bpAtPosition[j] = 0;
				seg->breakPointsi[seg->breakPoints[seg->numBp - 1]] = seg->breakPointsi[j];
				swap(seg->breakPoints, seg->breakPointsi[j], --seg->numBp);
				seg->breakPointsi[j] = -1;
			}
		}else if( seg->bpAtPosition[j] == 0 ){ // Else if elements are not adjacent and there is no breakpoint between them, place a breakpoint.
				seg->bpAtPosition[j] = 1;
				seg->breakPoints[seg->numBp] = j;
				seg->breakPointsi[j] = seg->numBp;
				seg->numBp ++;
		}

	}
	return seg;
}

int getLength(struct segment *seg){
	return seg->numBp;
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


