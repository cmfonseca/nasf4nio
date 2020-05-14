#include "problem.h"

struct problem *newProblem(const char *filename);
struct solution *setObjValue(double objv,struct solution *s);
int isSolEvaluated(struct solution *s);
