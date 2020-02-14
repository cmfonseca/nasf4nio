/* problem.h
 *
 * (C) 2019 Carlos M. Fonseca <cmfonsec@dei.uc.pt>
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

/* This header file contains all problem independent definitions */


  /*******************/
 /* Data structures */
/*******************/

struct problem;
struct solution;
struct move;
struct segment;


  /*************************/
 /* Problem instantiation */
/*************************/

/* newProblem() allocates a problem structure and initialises it. Being
 * problem-specific, the function arguments are deliberately left unspecified */
struct problem *newProblem();


  /*********************/
 /* Memory management */
/*********************/

struct solution *allocSolution(struct problem *p);
struct move *allocMove(struct problem *p);
struct segment *allocSegment(struct problem *p);

void freeProblem(struct problem *p);
void freeSolution(struct solution *s);
void freeMove(struct move *v);
void freeSegment(struct segment *g);


  /*******/
 /* I/O */
/*******/

void printProblem(struct problem *p);
void printSolution(struct solution *s);
void printMove(struct move *v);
void printSegment(struct segment *g);


  /***********************/
 /* Solution generation */
/***********************/

/* randomSolution() implements uniform random sampling of the solution space.
 * The input argument must be a pointer to a solution previously allocated with
 * allocSolution(), which is modified in place.
 * If randomMoveWOR() is implemented, randomSolution() must also initialise the
 * corresponding state by performing the equivalent to resetRandomMoveWOR().
 */
struct solution *randomSolution(struct solution *s);

/* randomNeighbour() implements uniform random sampling at a given distance
 * from a given solution, and modifies the solution in place. For d = 1, this
 * should be functionally equivalent to applyMove(s, randomMove(v, s))).
 * For d > 1, generating the solution directly may be more efficient than
 * applying a sequence of moves. Furthermore, generating random sequences of
 * moves that lead to uniform sampling at all distances may be simply imposible
 * depending on the neighbourhood structure.
 * If randomMoveWOR() is implemented, randomNeighbour() must also initialise the
 * corresponding state by performing the equivalent to resetRandomMoveWOR().
 */
struct solution *randomNeighbour(struct solution *s, int d);


  /***********************/
 /* Solution inspection */
/***********************/

/* getObjectiveValue() supports full and/or partial solution evaluation.
 * Once a solution is evaluated, results may be cached in the solution
 * itself so that it can be re-evaluated more efficiently after it is
 * modified by one or more calls to applyMove(). Therefore, the formal
 * argument is not const. Solution (re-)evaluation must occur before
 * this function returns, but the time at which actual evaluation occurs
 * is otherwise left unspecified.
 */
double getObjectiveValue(struct solution *s);

/* 
 * Number of direct neighbours of a given solution.
 */
int getNeighbourhoodSize(struct solution *s);

/* 
 * Distance to the farthest point in solution space.
 */
int getExcentricity(struct solution *s);


  /*******************/
 /* Move generation */
/*******************/

/* randomMove() implements uniform random sampling of the neighbourhood of a
 * given solution, with replacement.
 * The first input argument must be a pointer to a move previously allocated
 * with allocMove(), which is modified in place.
 */
struct move *randomMove(struct move *v, const struct solution *s);

/* randomMoveWOR() implements uniform random sampling of the neighbourhood of
 * a given solution, without replacement.
 * The first input argument must be a pointer to a move previously allocated
 * with allocMove(), which is modified in place. The function returns this
 * pointer if a new move is generated or NULL if there are no moves left.
 */
struct move *randomMoveWOR(struct move *v, struct solution *s);

/* resetRandomMoveWOR() resets the uniform random sampling without replacement
 * of the neighbourhood of a given solution, so that any move can be generated
 * by the next call to randomMoveWOR(). The function returns its input argument.
 */
struct solution *resetRandomMoveWOR(struct solution *s);


  /**************************/
 /* Operations on solutions*/
/**************************/

/* copySolution() copies the contents of the second argument to the first
 * argument, which must have been previously allocated with allocSolution().
 */
struct solution *copySolution(struct solution *dest, const struct solution *src);

/* applyMove() modifies a solution by applying a move to it. It is assumed
 * that the move was generated for, and possibly evaluated with respect   
 * to, that particular solution. The result of applying a move to a solution
 * other than that for which it was generated/evaluated (or a pristine copy of
 * it), including applying the same move to a solution more than once, is
 * undefined.
 * If randomMoveWOR() is implemented, applyMove must also reset the
 * corresponding state by performing the equivalent to resetRandomMoveWOR().
 */
struct solution *applyMove(struct solution *s, const struct move *v);


  /*******************/
 /* Move inspection */
/*******************/

/* getObjectiveIncrement() supports move evaluation with respect to the
 * solution for which it was generated, before it is actually applied to
 * that solution (if it ever is). The result of evaluating a move with
 * respect to a solution other than that for which it was generated 
 * (or to a pristine copy of it) is undefined.
 * Once a move is evaluated, results may be cached in the move itself, so
 * that they can be used by applyMove() to update the evaluation state of
 * the solution more efficiently.
 * In addition, results may also be cached in the solution in order
 * to speed up evaluation of future moves. Consequently, neither formal
 * argument is const.
 */
double getObjectiveIncrement(struct move *v, struct solution *s);


  /**********************/
 /* Segment operations */
/**********************/

struct segment *initSegment(struct segment *g, const struct solution *s1, const struct solution *s2);

/* randomSolutionInSegment() implements uniform random sampling of the
 * solutions in a segment, including the extremes, with replacement (aka
 * uniform crossover).
 * The first input argument must have been previously allocated with
 * allocSolution(), and is modified in place.
 * TBC: Can partial evaluation information be propagated to the new solution?
 */
struct solution *randomSolutionInSegment(struct solution *s, const struct segment *g);

/* randomNeighbourTowards() samples the solutions in a segment that are at a
 * given distance from the first extreme uniformly at random, with replacement.
 * The input argument must be a pointer to a solution previously allocated with
 * allocSolution(), which is modified in place.
 * TBC: Can partial evaluation information be propagated to the new solution?
 */
struct solution *randomNeighbourTowards(struct solution *s, const struct segment *g, int d);

/* randomMoveTowards() implements uniform random sampling with replacement
 * of the part of the neighbourhood of the first extreme of a segment that
 * lies in the segment itself.
 * The first input argument must be a pointer to a move previously allocated
 * with allocMove(), which is modified in place. The segment is not modified.
 * FIXME: should struct segment be const?
 */
struct move *randomMoveTowards(struct move *v, struct segment *g);

struct move *randomMoveTowardsWOR(struct move *v, struct segment *g);

/* randomMoveAway() implements uniform random sampling with replacement
 * of the part of the neighbourhood of the first extreme of a segment that
 * lies at a greater distance from the second extreme than the first extreme
 * itself.
 * The first input argument must be a pointer to a move previously allocated
 * with allocMove(), which is modified in place. The segment is not modified.
 */
struct move *randomMoveAway(struct move *v, struct segment *g);

struct move *randomMoveAwayWOR(struct move *v, struct segment *g);

#if 0

/* It is unclear whether there will be a use case for this, but leave it
 * commented out here for completeness.
 */

/* randomMoveAside() implements uniform random sampling with replacement
 * of the part of the neighbourhood of the first extreme of a segment that
 * lies at exactly the same distance from the second extreme than the first
 * extreme itself. This may be an empty set for some points or neighbourhoods,
 * in which case the function returns NULL.
 * The first input argument must be a pointer to a move previously allocated
 * with allocMove(), which is modified in place. The segment is not modified.
 */
struct move *randomMoveAside(struct move *v, const struct segment *g);

struct move *randomMoveAsideWOR(struct move *v, struct segment *g);

#endif

struct segment *applyMoveToSegment(struct segment *g, const struct move *v);

/* returns the length of the segment. The segment is not const so that
 * the value and/or any other intermediate results can be cached in the
 * segment if desired.
 */
int getLength(struct segment *g);



