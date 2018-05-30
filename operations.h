#ifndef OPERATIONS
#define OPERATIONS

#include "solver.h"

double e2_scalar_product_array (int n, const double *u, const double *v);
double e2_norm_array (int n, const double *u);
void copy_array (int n, const double *src, double *dst);
void multiply_matrix_vector (int n, msr_matrix_rhs *matrix_rhs_storage, double *x, double *y);

void linear_combination_1 (int n, const double *x, const double *y, const double alpha,
                           double *res);

void linear_combination_5 (int n,    // array size
                      const double *x,  // input array X
                      const double *y,  // input array Y
                      double alpha,     // multiplier
                      double beta,      // multiplier
                      double *res       // input/output array
                      );

void linear_combination_2 (int n, const double *x, const double *y,
                           double alpha, double beta,
                           double *res);


void jacobi_preconditioner (int n, msr_matrix_rhs *matrix_rhs_storage, double *r, double *z);



#endif
