#include "operations.h"
#include "solver.h"
#include <math.h>

double e2_scalar_product_array (int n, const double *u, const double *v)
{
  double sum = 0.;

  for (int i = 0; i < n; ++i)
    {
      sum += u[i] * v[i];
    }

  return sum;
}

double e2_norm_array (int n, const double *u)
{
  return sqrt (e2_scalar_product_array (n, u, u));
}

void copy_array (int n, const double *src, double *dst)
{
  for (int i = 0; i < n; ++i)
    dst[i] = src[i];
}

void jacobi_preconditioner (int n, msr_matrix_rhs *matrix_rhs_storage, double *r, double *z)
{
  for (int i = 0; i < n; ++i)
    z[i] = r[i] / matrix_rhs_storage->matrix[i];
}

void multiply_matrix_vector (int n, msr_matrix_rhs *matrix_rhs_storage, double *x, double *y)
{
  const double *a = matrix_rhs_storage->matrix;
  const int *indexes = matrix_rhs_storage->row_non_zeros;
  int i, j, s, len;
  double sum;

  for (i = 0; i < n; i++)
    {
      sum = a[i] * x[i];
      len = indexes[i + 1] - indexes[i];
      s = indexes[i];

      for (j = 0; j < len; j++)
        sum += a[s + j] * x[indexes[s + j]];

      y[i] = sum;
    }
}

void linear_combination_1 (int n, const double *x, const double *y, const double alpha,
                           double *res)
{
  for (int i = 0; i < n; ++i)
    {
      res[i] = x[i] + alpha * y[i];
    }
}

void linear_combination_5 (int n,    // array size
                      const double *x,  // input array X
                      const double *y,  // input array Y
                      double alpha,     // multiplier
                      double beta,      // multiplier
                      double *res       // input/output array
                      )
{
  for (int i = 0; i < n; ++i)
    {
      res[i] = alpha * x[i] + beta * y[i];
    }
}

void linear_combination_2 (int n, const double *x, const double *y,
                           double alpha, double beta,
                           double *res)
{
  for (int i = 0; i < n; ++i)
    {
      res[i] = x[i] + beta * (res[i] + alpha * y[i]);
    }
}
