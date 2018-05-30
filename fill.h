#ifndef FILL_H
#define FILL_H


#include "solver.h"
#include "parser.h"

int fill_by_prev_layer (solver_t *solver);

int fill_dot (solver_t *solver, int dot_number);

int fill_dot_0 (solver_t *solver, int dot_number);
int fill_dot_1 (solver_t *solver, int dot_number);
int fill_dot_2 (solver_t *solver, int dot_number);
int fill_dot_3 (solver_t *solver, int dot_number);
int fill_dot_4 (solver_t *solver, int dot_number);
int fill_dot_5 (solver_t *solver, int dot_number);

typedef struct
{
  double value;
  int column;
} msr_pair;

double f0 ();
double f1 ();
double f2 ();


#endif
