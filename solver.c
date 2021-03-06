#include "solver.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "it_algoritm.h"
#include "neighbors_and_type_init.h"


void solver_init (solver_t *solver,
                 solver_ws *base_ws_arg,
                 msr_matrix_rhs *matrix_rhs_storage_arg,
                 msr_matrix_rhs *matrix_rhs_storage_preconditioned_arg,
                 int MX_arg,
                 int MY_arg,
                 int N_arg,
                 double T_arg,
                 double border_omega_arg,
                 double mu_arg,
                 int solver_type_arg)
{
  int total_dots_amount = (MY_arg * 3 + 1) * (MX_arg + 1) + MX_arg * (MY_arg + 1);
  int total_equations_amount = 3 * total_dots_amount;

  solver->MX = MX_arg;
  solver->MY = MY_arg;
  solver->N = N_arg;
  solver->T = T_arg;
  solver->hx = M_PI / solver->MX;
  solver->hy = M_PI / solver->MY;
  solver->tau = solver->T / solver->N;

  solver->border_omega = border_omega_arg;
  solver->mu = mu_arg;

  solver->solver_type = solver_type_arg;

  solver->base_ws = base_ws_arg;

  solver->base_ws->prev_answer = (double *)malloc(total_equations_amount * sizeof (double));

  solver->base_ws->total_dots_amount = total_dots_amount;
  solver->base_ws->total_equation_amount = total_equations_amount;

  solver->base_ws->matrix_rhs_storage = matrix_rhs_storage_arg;
  solver->base_ws->matrix_rhs_storage_preconditioned = matrix_rhs_storage_preconditioned_arg;

  msr_matrix_rhs_init (total_equations_amount, solver->base_ws->matrix_rhs_storage);
  msr_matrix_rhs_init (total_equations_amount, solver->base_ws->matrix_rhs_storage_preconditioned);

  solver->base_ws->workspace = (double *)malloc(total_equations_amount * sizeof (double) * 7);

  if (solver->solver_type == 1)
    {
      // Laspack init
    }
  else
    {
      // Own ws init
    }

  init_neighbors_and_type (solver);

}

void msr_matrix_rhs_init (int total_equations_amount, msr_matrix_rhs *matrix_rhs_storage)
{
  matrix_rhs_storage->n = total_equations_amount;
  matrix_rhs_storage->matrix = (double *)malloc((total_equations_amount + total_equations_amount * 9 + 1) * sizeof (double));
  matrix_rhs_storage->rhs = (double *)malloc(total_equations_amount * sizeof (double));
  matrix_rhs_storage->row_non_zeros = (int *)malloc((total_equations_amount + total_equations_amount * 9 + 1 + 1) * sizeof (int));

  memset (matrix_rhs_storage->row_non_zeros, 0, (total_equations_amount + total_equations_amount * 9 + 1) * sizeof (int));
  matrix_rhs_storage->row_non_zeros[0] = matrix_rhs_storage->n + 1;
}

void solver_free (solver_t *solver)
{
  free (solver->top_neighbor);
  free (solver->bottom_neighbor);
  free (solver->left_neighbor);
  free (solver->right_neighbor);
  free (solver->dot_type);

  free (solver->base_ws->matrix_rhs_storage->matrix);
  free (solver->base_ws->matrix_rhs_storage->rhs);
  free (solver->base_ws->matrix_rhs_storage->row_non_zeros);

  free (solver->base_ws->matrix_rhs_storage_preconditioned->matrix);
  free (solver->base_ws->matrix_rhs_storage_preconditioned->rhs);
  free (solver->base_ws->matrix_rhs_storage_preconditioned->row_non_zeros);

  free (solver->base_ws->workspace);
}


int get_equation_g (int dot_number)
{
  return dot_number * 3;
}

int get_equation_v1 (int dot_number)
{
  return dot_number * 3 + 1;
}

int get_equation_v2 (int dot_number)
{
  return dot_number * 3 + 2;
}

int solver_run (solver_t *solver)
{
  int iters = bcgstab (solver);
  printf ("Iters = %d\n", iters);
  return iters;
}
