/* solver.c
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

#include "problem.h"

/* This header file contains all solver independent definitions */

/* Data structure */
struct solverState;

/* Solver instantiation - function arguments are purposely left unspecified,
 * but they should include a problem instance
 */
struct solverState *newSolver();

/* Memory management */
void freeSolver(struct solverState *ss);

/* I/O */
void printSolverState(struct solverState *ss);

/* State iterator */
struct solverState *nextSolverState(struct solverState *ss);

/* State inspection */
struct solution *getSolverSolution(struct solverState *ss);

