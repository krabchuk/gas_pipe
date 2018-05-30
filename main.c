#include "solver.h"
#include "it_algoritm.h"
#include "fill.h"
#include <stdio.h>

int main(void)
{
  solver_t solver;
  solver_ws base_ws;
  msr_matrix_rhs matrix_rhs_storage;
  msr_matrix_rhs matrix_rhs_storage_preconditioned;

  solver_init (&solver,
               &base_ws,
               &matrix_rhs_storage,
               &matrix_rhs_storage_preconditioned,
               3,
               3,
               10,
               1,
               1,
               1,
               1);

  fill_zero_layer (&solver);
  fill_by_prev_layer (&solver);

  for (int i = 0; i < solver.base_ws->matrix_rhs_storage->total_nz; ++i)
    printf ("%d\n", solver.base_ws->matrix_rhs_storage->row_non_zeros[i]);
#if 0
  {1, 1, 1, 0,
   6, 2, 3, 0,
   0, 2, 5, 8,
   0, 0, 4, 6};
#endif
  solver.base_ws->matrix_rhs_storage->n = 4;
  solver.base_ws->matrix_rhs_storage->rhs[0] = 1;
  solver.base_ws->matrix_rhs_storage->rhs[1] = 2;
  solver.base_ws->matrix_rhs_storage->rhs[2] = 5;
  solver.base_ws->matrix_rhs_storage->rhs[3] = 6;

  solver.base_ws->matrix_rhs_storage->matrix[0] = 1;
  solver.base_ws->matrix_rhs_storage->matrix[1] = 2;
  solver.base_ws->matrix_rhs_storage->matrix[2] = 5;
  solver.base_ws->matrix_rhs_storage->matrix[3] = 6;


  solver.base_ws->matrix_rhs_storage->row_non_zeros[0] = 5;
  solver.base_ws->matrix_rhs_storage->row_non_zeros[1] = 5;
  solver.base_ws->matrix_rhs_storage->row_non_zeros[2] = 5;
  solver.base_ws->matrix_rhs_storage->row_non_zeros[3] = 5;

  if (solver_run (&solver) < 0)
    return 0;

  for (int i = 0; i < 4; ++i)
    printf ("%e\n", solver.base_ws->prev_answer[i]);

//  solver_free (&solver);

  return 0;
}
