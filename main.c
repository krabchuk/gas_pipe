#include "solver.h"

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
               2,
               1,
               1,
               1,
               1,
               1);

  solver_free (&solver);

  return 0;
}
