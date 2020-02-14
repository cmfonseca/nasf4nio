/* sga.h
 *
 * (C) 2018 Carlos M. Fonseca <cmfonsec@dei.uc.pt>
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

#include "solver.h"

struct solverState *newSolver(struct problem *p, int popsize);

/* Solver specific functions */

double getSelectivePressure(struct solverState *ss);
double getRecombinationRate(struct solverState *ss);
double getMutationRate(struct solverState *ss);

void setSelectivePressure(struct solverState *ss, double SP);
void setRecombinationRate(struct solverState *ss, double Px);
void setMutationRate(struct solverState *ss, double Pm);

