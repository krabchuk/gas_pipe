#include "it_algoritm.h"
#include "solver.h"
#include "operations.h"
#include <stdio.h>
#include <math.h>


int bcgstab (solver_t *solver)
{
  int n = solver->base_ws->matrix_rhs_storage->n;
  double *rhs = solver->base_ws->matrix_rhs_storage->rhs;
  unsigned int itmax = 1000;
  double itprec = 1e-9;
  msr_matrix_rhs *matrix = solver->base_ws->matrix_rhs_storage;
  double *u = solver->base_ws->prev_answer;
  double rhs_norm;              // initial norm of right hand side
  double residual_norm;         // norm of current residual
  double stop_test;             // stopping criterion
  unsigned int iter;            // Output: number of iterations proceeded
  double masheps = 1.e-16;      // ANSI/IEEE 754-1985 floating point precision
  double min_for_division = 1.e-64;     // Minimal for division

  double c_1, alpha, beta, rho, omega;  // Algorithm variables
  double scalp;                 // Temporary variable

  // initialize workspace arrays
  double *r = solver->base_ws->workspace;
  double *p = r + n;
  double *s = p + n;
  double *t = s + n;
  double *v = t + n;
  double *w = v + n;
  double *z = w + n;

  // Zero initial approximation
  // r = rhs
  copy_array (n, rhs, r);
  // Compute norm of residual
  residual_norm = e2_norm_array (n, r);
  // Compute norm of right hand side
  rhs_norm = residual_norm;

  if (rhs_norm < masheps)
    rhs_norm = masheps;

  stop_test = residual_norm / rhs_norm;
  if (stop_test < itprec)
    return 0;

#if PRINT
  printf ("\tit=%3d, residual = %le, rhs = %le\n", 0, residual_norm, rhs_norm);
#endif

  iter = 1;
  jacobi_preconditioner (n, matrix, r, z);


  // Initialize r: r = z
  copy_array (n, z, r);

  // Initialize p: p = z
  copy_array (n, z, p);

  // Compute rho = (z, r)
  rho = e2_scalar_product_array (n, z, r);
  if (fabs (rho) < min_for_division)
    {
      // Method cannot be applied with this initial guesss
      printf ("Fatal error: acceleration algorithm breaks down\n");
      return -1;
    }

  // Iterative loop
  for (; iter <= itmax; iter++)
    {
      // Compute v = B^{-1}A(p)
      // First w = A(p)
      multiply_matrix_vector (n, matrix, p, w);
      // Next v = B^{-1} w
      jacobi_preconditioner (n, matrix, w, v);

      // Compute alpha = rho / (z, v)
      // Compute first c_1 = (z, v)
      c_1 = e2_scalar_product_array (n, z, v);
      if (fabs (c_1) < min_for_division)
        {
          // Method cannot be applied with this initial guesss
          printf ("Fatal error: acceleration algorithm breaks down\n");
          return -2;
        }
      // Compute alpha = rho / (z, v)
      alpha = rho / c_1;

      // s = r - alpha v
      linear_combination_1 (n, r, v, -alpha, s);

      // Compute t = B^{-1}A(s)
      // First w = A(s)
      multiply_matrix_vector (n, matrix, s, w);

      // Next t = B^{-1} w
      jacobi_preconditioner (n, matrix, w, t);

      // Compute omega = (t, s) / (t, t)
      scalp = e2_scalar_product_array (n, t, t);
      if (fabs (scalp) < min_for_division)
        {
          // t = B^{-1}A(s) = 0 => s = 0 that is r = alpha v = alpha B^{-1}A(p)
          // so u = u + alpha is precise solution
          // u = u + alpha p + omega s
          linear_combination_5 (n, p, s, alpha, 0., u);
          break;
        }
      // Compute omega = (t, s) / (t, t)
      omega = e2_scalar_product_array (n, t, s) / scalp;
      if (fabs (omega) < min_for_division)
        {
          // Method cannot be applied with this initial guesss
          printf ("Fatal error: acceleration algorithm breaks down\n");
          return -4;
        }

      // u = u + alpha p + omega s
      linear_combination_5 (n, p, s, alpha, omega, u);

      // r = s - omega t
      linear_combination_1 (n, s, t, -omega, r);

      // Compute residual norm
      residual_norm = e2_norm_array (n, r);

      // Check stopping criterion
      stop_test = residual_norm / rhs_norm;
#if PRINT
      printf ("\tit=%3d, residual = %le\n", iter, residual_norm);
#endif
      if (stop_test < itprec)
        break;

      // Compute rho_next = c_1 = (z, r)
      c_1 = e2_scalar_product_array (n, z, r);
      if (fabs (c_1) < min_for_division)
        {
          // Method cannot be applied with this initial guesss
          printf ("Fatal error: acceleration algorithm breaks down\n");
          return -5;
        }
      // beta = (rho_next / \rho) / (alpha / beta)
      beta = c_1 / rho * alpha / omega;
      // rho = rho_next
      rho = c_1;

      // p = r + beta (p - omega v)
      linear_combination_2 (n, r, v, - omega, beta, p);
    }

  if (iter >= itmax)
    {
      // Convergence too slow or failed
      printf ("Fatal error: failure to converge in max %d iterations\n", itmax);

      return -10;
    }

  return iter;
}


