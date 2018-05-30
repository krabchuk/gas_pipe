#include "fill.h"
#include <stdio.h>
#include <math.h>

double f0 ()
{
  return 0;
}

double f1 ()
{
  return 0;
}

double f2 ()
{
  return 0;
}

int fill_by_prev_layer (solver_t *solver)
{
  int res = 0;

  for (int i = 0; i < solver->base_ws->total_dots_amount; ++i)
    {
      res = fill_dot (solver, i);
      if (res)
        {
          printf ("Error %d in dot %d\n", res, i);
          return -1;
        }
    }
  return 0;
}

int fill_dot (solver_t *solver, int dot_number)
{
  int res = 0;
  switch (solver->dot_type[dot_number])
    {
    case 0:
      res = fill_dot_0 (solver, dot_number);
      break;
    case 1:
      res = fill_dot_1 (solver, dot_number);
      break;
    case 2:
      res = fill_dot_2 (solver, dot_number);
      break;
    case 3:
      res = fill_dot_3 (solver, dot_number);
      break;
    case 4:
      res = fill_dot_4 (solver, dot_number);
      break;
    case 5:
      res = fill_dot_5 (solver, dot_number);
      break;
    default:
      printf ("Unknown dot type\n");
      res = -1;
    }

  return res;
}

int fill_dot_0 (solver_t *solver, int dot_number)
{
  int top = solver->top_neighbor[dot_number];
  int bot = solver->bottom_neighbor[dot_number];
  int left = solver->left_neighbor[dot_number];
  int right = solver->right_neighbor[dot_number];
  
  int eq_g = get_equation_g (dot_number);
  int eq_v1 = get_equation_v1 (dot_number);
  int eq_v2 = get_equation_v2 (dot_number);
  
  double *prev = solver->base_ws->prev_answer;
  
  double *matrix = solver->base_ws->matrix_rhs_storage->matrix;
  double *rhs = solver->base_ws->matrix_rhs_storage->rhs;
  int *row_non_zeros = solver->base_ws->matrix_rhs_storage->row_non_zeros;
  int n = solver->base_ws->matrix_rhs_storage->n;
  
  double tau = solver->tau;
  double hx = solver->hx;
  double hy = solver->hy;

  double mu = solver->mu;

  double exp_minusG = exp (-prev[eq_g]);

  msr_pair pairs[8];

  double mu_voln = prev[eq_g];

  for (int i = 0; i < n; i += 3)
    if (prev[i] < mu_voln)
      mu_voln = prev[i];

  mu_voln = mu * exp (-mu_voln);
  

  
  
  // Equation for G
  solver->base_ws->matrix_rhs_storage->total_nz += 9;
  
  rhs[eq_g] = 4 * prev[eq_g] + tau * prev[eq_g] * 
      ((prev[get_equation_v1 (right)] - prev[get_equation_v1 (left)]) / hx +
      (prev[get_equation_v2 (top)] - prev[get_equation_v2 (bot)]) / hy) + 4 * tau * f0 ();
  
  row_non_zeros[eq_g + 1] = row_non_zeros[eq_g] + 8;
  
  matrix[eq_g] = 4;
  
  pairs[0].value = - (prev[eq_v1] + prev[get_equation_v1 (left)]) * tau / hx;
  pairs[0].column = get_equation_g (left);
  
  pairs[1].value = (prev[eq_v1] + prev[get_equation_v1 (right)]) * tau / hx;
  pairs[1].column = get_equation_g (right);
  
  pairs[2].value = - (prev[eq_v2] + prev[get_equation_v2 (bot)]) * tau / hy;
  pairs[2].column = get_equation_g (bot);
  
  pairs[3].value = (prev[eq_v2] + prev[get_equation_v2 (top)]) * tau / hy;
  pairs[3].column = get_equation_g (top);
  
  pairs[4].value = -2 * tau / hx;
  pairs[4].column = get_equation_v1 (left);
  
  pairs[5].value = 2 * tau / hx;
  pairs[5].column = get_equation_v1 (right);
  
  pairs[6].value = -2 * tau / hy;
  pairs[6].column = get_equation_v2 (bot);
  
  pairs[7].value = 2 * tau / hy;
  pairs[7].column = get_equation_v2 (top);
  
  for (int i = 0; i < 8; ++i)
    {
      matrix[row_non_zeros[eq_g] + i] = pairs[i].value;
      row_non_zeros[row_non_zeros[eq_g] + i] = pairs[i].column;
    }
  
  // Ecuation for V1
  solver->base_ws->matrix_rhs_storage->row_non_zeros += 7;

  rhs[eq_v1] = 6 * prev[eq_v1] + 6 * tau * f1 () + 1.5 * tau * prev[eq_v1] * (prev[get_equation_v2 (top)] - prev[get_equation_v2 (bot)]) / hy -
      (mu_voln - mu * exp_minusG) * 6 * tau * 4 * (prev[get_equation_v1 (right)] - 2 * prev[eq_v1] + prev[get_equation_v1 (left)]) / (3 * hx * hx) +
      (mu_voln - mu * exp_minusG) * 6 * tau * (prev[get_equation_v1 (top)] - 2 * prev[eq_v1] + prev[get_equation_v1 (bot)]) / (hy * hy) +
      mu * exp_minusG * (prev[get_equation_v2 (solver->top_neighbor[right])] - prev[get_equation_v2 (solver->bottom_neighbor[right])]
      - prev[get_equation_v2 (solver->top_neighbor[left])] + prev[get_equation_v2 (solver->bottom_neighbor[left])]);

  row_non_zeros[eq_v1 + 1] = row_non_zeros[eq_v1] + 6;

  matrix[eq_v1] = 6 + 4 * tau * mu_voln * (4 / (hx * hx) + 3 / (hy * hy));

  pairs[0].value = -tau * (prev[eq_v1] + prev[get_equation_v1 (left)]) / hx - mu_voln * 8 * tau / (hx * hx);
  pairs[0].column = get_equation_v1 (left);

  pairs[1].value = tau * (prev[eq_v1] + prev[get_equation_v1 (right)]) / hx - mu_voln * 8 * tau / (hx * hx);
  pairs[1].column = get_equation_v1 (right);

  pairs[2].value = -1.5 * tau * (prev[eq_v2] + prev[get_equation_v2 (bot)]) / hy - mu_voln * 6 * tau / (hy * hy);
  pairs[2].column = get_equation_v1 (bot);

  pairs[3].value = 1.5 * tau * (prev[eq_v2] + prev[get_equation_v2 (top)]) / hy - mu_voln * 6 * tau / (hy * hy);
  pairs[3].column = get_equation_v1 (top);

  pairs[4].value = -3 * tau / hx;
  pairs[4].column = get_equation_g (left);

  pairs[5].value = 3 * tau / hx;
  pairs[5].value = get_equation_g (right);

  for (int i = 0; i < 6; ++i)
    {
      matrix[row_non_zeros[eq_v1] + i] = pairs[i].value;
      row_non_zeros[row_non_zeros[eq_v1] + i] = pairs[i].column;
    }

  // Equation for V2
  solver->base_ws->matrix_rhs_storage->row_non_zeros += 7;

  rhs[eq_v2] = 6 * prev[eq_v2] + 6 * tau * f2 () + 1.5 * tau * prev[eq_v2] * (prev[get_equation_v1 (right)] - prev[get_equation_v1 (left)]) / hx -
      (mu_voln - mu * exp_minusG) * 6 * tau * 4 * (prev[get_equation_v2 (top)] - 2 * prev[eq_v2] + prev[get_equation_v2 (bot)]) / (3 * hy * hy) +
      (mu_voln - mu * exp_minusG) * 6 * tau * (prev[get_equation_v2 (right)] - 2 * prev[eq_v2] + prev[get_equation_v2 (left)]) / (hx * hx) +
      mu * exp_minusG * (prev[get_equation_v1 (solver->bottom_neighbor[left])] - prev[get_equation_v1 (solver->top_neighbor[left])]
      - prev[get_equation_v1 (solver->bottom_neighbor[right])] + prev[get_equation_v2 (solver->top_neighbor[right])]);

  row_non_zeros[eq_v2 + 1] = row_non_zeros[eq_v2] + 6;

  matrix[eq_v2] = 6 + 4 * tau * mu_voln * (4 / (hy * hy) + 3 / (hx * hx));

  pairs[0].value = -tau * (prev[eq_v2] + prev[get_equation_v2 (bot)]) / hy - mu_voln * 8 * tau / (hy * hy);
  pairs[0].column = get_equation_v2 (bot);

  pairs[1].value = tau * (prev[eq_v2] + prev[get_equation_v2 (top)]) / hy - mu_voln * 8 * tau / (hy * hy);
  pairs[1].column = get_equation_v2 (top);

  pairs[2].value = -1.5 * tau * (prev[eq_v1] + prev[get_equation_v1 (left)]) / hx - mu_voln * 6 * tau / (hx * hx);
  pairs[2].column = get_equation_v2 (left);

  pairs[3].value = 1.5 * tau * (prev[eq_v1] + prev[get_equation_v1 (right)]) / hx - mu_voln * 6 * tau / (hx * hx);
  pairs[3].column = get_equation_v2 (right);

  pairs[4].value = -3 * tau / hy;
  pairs[4].column = get_equation_g (bot);

  pairs[5].value = 3 * tau / hy;
  pairs[5].value = get_equation_g (top);

  for (int i = 0; i < 6; ++i)
    {
      matrix[row_non_zeros[eq_v2] + i] = pairs[i].value;
      row_non_zeros[row_non_zeros[eq_v2] + i] = pairs[i].column;
    }


  return 0;
}


int fill_dot_1 (solver_t *solver, int dot_number)
{
//  int top = solver->top_neighbor[dot_number];
//  int bot = solver->bottom_neighbor[dot_number];
//  int left = solver->left_neighbor[dot_number];
  int right = solver->right_neighbor[dot_number];

  int eq_g = get_equation_g (dot_number);
  int eq_v1 = get_equation_v1 (dot_number);
  int eq_v2 = get_equation_v2 (dot_number);

  double *prev = solver->base_ws->prev_answer;

  double *matrix = solver->base_ws->matrix_rhs_storage->matrix;
  double *rhs = solver->base_ws->matrix_rhs_storage->rhs;
  int *row_non_zeros = solver->base_ws->matrix_rhs_storage->row_non_zeros;
//  int n = solver->base_ws->matrix_rhs_storage->n;

  double tau = solver->tau;
  double hx = solver->hx;
//  double hy = solver->hy;

//  double mu = solver->mu;

//  double exp_minusG = exp (-prev[eq_g]);

  msr_pair pairs[8];

//  double mu_voln = prev[eq_g];

//  for (int i = 0; i < n; i += 3)
//    if (prev[i] < mu_voln)
//      mu_voln = prev[i];

//  mu_voln = mu * exp (-mu_voln);

  // Equation for G
   solver->base_ws->matrix_rhs_storage->total_nz += 4;

   rhs[eq_g] = 2 * prev[eq_g] + tau * (prev[get_equation_v1 (right)] - prev[eq_v1]) / hx + 2 * tau * f0 () +
       tau * (prev[eq_g] * prev[eq_v1] - 2.5 * prev[get_equation_g (right)] * prev[get_equation_v1 (right)]
       + 2 * prev[get_equation_g (solver->right_neighbor[right])] * prev[get_equation_v1 (solver->right_neighbor[right])] -
       0.5 * prev[get_equation_g (solver->right_neighbor[solver->right_neighbor[right]])] * prev[get_equation_v1 (solver->right_neighbor[solver->right_neighbor[right]])] +
       (2 - prev[eq_g]) * (prev[eq_v1] - 2.5 * prev[get_equation_v1 (right)] + 2 * prev[get_equation_v1 (solver->right_neighbor[right])] - 0.5 * prev[get_equation_v1 (solver->right_neighbor[solver->right_neighbor[right]])])) / hx;

   row_non_zeros[eq_g + 1] = row_non_zeros[eq_g] + 3;

   matrix[eq_g] = 2 - tau * prev[eq_v1] / hx;

   pairs[0].value = tau * prev[get_equation_v1 (right)] / hx;
   pairs[0].column = get_equation_g (right);

   pairs[1].value = 2 * tau / hx;
   pairs[1].column = get_equation_v1 (right);

   pairs[2].value = -2 * tau / hx;
   pairs[2].column = eq_v1;

   for (int i = 0; i < 3; ++i)
     {
       matrix[row_non_zeros[eq_g] + i] = pairs[i].value;
       row_non_zeros[row_non_zeros[eq_g] + i] = pairs[i].column;
     }

   //Equation for V1
   solver->base_ws->matrix_rhs_storage->total_nz += 1;

   row_non_zeros[eq_v1 + 1] = row_non_zeros[eq_v1];

   rhs[eq_v1] = solver->border_omega;

   matrix[eq_v1] = 1;

   //Equation for V2
   solver->base_ws->matrix_rhs_storage->total_nz += 1;

   row_non_zeros[eq_v2 + 1] = row_non_zeros[eq_v2];

   rhs[eq_v2] = 0;

   matrix[eq_v2] = 1;

   return 0;
}

int fill_dot_4 (solver_t *solver, int dot_number)
{
//  int top = solver->top_neighbor[dot_number];
//  int bot = solver->bottom_neighbor[dot_number];
//  int left = solver->left_neighbor[dot_number];
  int right = solver->right_neighbor[dot_number];

  int eq_g = get_equation_g (dot_number);
  int eq_v1 = get_equation_v1 (dot_number);
  int eq_v2 = get_equation_v2 (dot_number);

  double *prev = solver->base_ws->prev_answer;

  double *matrix = solver->base_ws->matrix_rhs_storage->matrix;
  double *rhs = solver->base_ws->matrix_rhs_storage->rhs;
  int *row_non_zeros = solver->base_ws->matrix_rhs_storage->row_non_zeros;
//  int n = solver->base_ws->matrix_rhs_storage->n;

  double tau = solver->tau;
  double hx = solver->hx;
//  double hy = solver->hy;

//  double mu = solver->mu;

//  double exp_minusG = exp (-prev[eq_g]);

  msr_pair pairs[8];

//  double mu_voln = prev[eq_g];

//  for (int i = 0; i < n; i += 3)
//    if (prev[i] < mu_voln)
//      mu_voln = prev[i];

//  mu_voln = mu * exp (-mu_voln);

  // Equation for G
   solver->base_ws->matrix_rhs_storage->total_nz += 4;

   rhs[eq_g] = 2 * prev[eq_g] + tau * (prev[get_equation_v1 (right)] - prev[eq_v1]) / hx + 2 * tau * f0 () +
       tau * (prev[eq_g] * prev[eq_v1] - 2.5 * prev[get_equation_g (right)] * prev[get_equation_v1 (right)]
       + 2 * prev[get_equation_g (solver->right_neighbor[right])] * prev[get_equation_v1 (solver->right_neighbor[right])] -
       0.5 * prev[get_equation_g (solver->right_neighbor[solver->right_neighbor[right]])] * prev[get_equation_v1 (solver->right_neighbor[solver->right_neighbor[right]])] +
       (2 - prev[eq_g]) * (prev[eq_v1] - 2.5 * prev[get_equation_v1 (right)] + 2 * prev[get_equation_v1 (solver->right_neighbor[right])] - 0.5 * prev[get_equation_v1 (solver->right_neighbor[solver->right_neighbor[right]])])) / hx;

   row_non_zeros[eq_g + 1] = row_non_zeros[eq_g] + 3;

   matrix[eq_g] = 2 - tau * prev[eq_v1] / hx;

   pairs[0].value = tau * prev[get_equation_v1 (right)] / hx;
   pairs[0].column = get_equation_g (right);

   pairs[1].value = 2 * tau / hx;
   pairs[1].column = get_equation_v1 (right);

   pairs[2].value = -2 * tau / hx;
   pairs[2].column = eq_v1;

   for (int i = 0; i < 3; ++i)
     {
       matrix[row_non_zeros[eq_g] + i] = pairs[i].value;
       row_non_zeros[row_non_zeros[eq_g] + i] = pairs[i].column;
     }

   //Equation for V1
   solver->base_ws->matrix_rhs_storage->total_nz += 1;

   row_non_zeros[eq_v1 + 1] = row_non_zeros[eq_v1];

   rhs[eq_v1] = 0;

   matrix[eq_v1] = 1;

   //Equation for V2
   solver->base_ws->matrix_rhs_storage->total_nz += 1;

   row_non_zeros[eq_v2 + 1] = row_non_zeros[eq_v2];

   rhs[eq_v2] = 0;

   matrix[eq_v2] = 1;

   return 0;
}

int fill_dot_5 (solver_t *solver, int dot_number)
{
//  int top = solver->top_neighbor[dot_number];
//  int bot = solver->bottom_neighbor[dot_number];
  int left = solver->left_neighbor[dot_number];
//  int right = solver->right_neighbor[dot_number];

  int eq_g = get_equation_g (dot_number);
  int eq_v1 = get_equation_v1 (dot_number);
  int eq_v2 = get_equation_v2 (dot_number);

  double *prev = solver->base_ws->prev_answer;

  double *matrix = solver->base_ws->matrix_rhs_storage->matrix;
  double *rhs = solver->base_ws->matrix_rhs_storage->rhs;
  int *row_non_zeros = solver->base_ws->matrix_rhs_storage->row_non_zeros;
//  int n = solver->base_ws->matrix_rhs_storage->n;

  double tau = solver->tau;
  double hx = solver->hx;
//  double hy = solver->hy;

//  double mu = solver->mu;

//  double exp_minusG = exp (-prev[eq_g]);

  msr_pair pairs[8];

//  double mu_voln = prev[eq_g];

//  for (int i = 0; i < n; i += 3)
//    if (prev[i] < mu_voln)
//      mu_voln = prev[i];

//  mu_voln = mu * exp (-mu_voln);

  // Equation for G
   solver->base_ws->matrix_rhs_storage->total_nz += 4;

   rhs[eq_g] = 2 * prev[eq_g] + tau * prev[eq_g] * (prev[eq_v1] - prev[get_equation_v1 (left)]) / hx + 2 * tau * f0 () -
       tau * (prev[eq_g] * prev[eq_v1] - 2.5 * prev[get_equation_g (left)] * prev[get_equation_v1 (left)]
       + 2 * prev[get_equation_g (solver->left_neighbor[left])] * prev[get_equation_v1 (solver->left_neighbor[left])] -
       0.5 * prev[get_equation_g (solver->left_neighbor[solver->left_neighbor[left]])] * prev[get_equation_v1 (solver->left_neighbor[solver->left_neighbor[left]])] +
       (2 - prev[eq_g]) * (prev[eq_v1] - 2.5 * prev[get_equation_v1 (left)] + 2 * prev[get_equation_v1 (solver->left_neighbor[left])] - 0.5 * prev[get_equation_v1 (solver->left_neighbor[solver->left_neighbor[left]])])) / hx;

   row_non_zeros[eq_g + 1] = row_non_zeros[eq_g] + 3;

   matrix[eq_g] = 2 + tau * prev[eq_v1] / hx;

   pairs[0].value = -tau * prev[get_equation_v1 (left)] / hx;
   pairs[0].column = get_equation_g (left);

   pairs[1].value = -2 * tau / hx;
   pairs[1].column = get_equation_v1 (left);

   pairs[2].value = 2 * tau / hx;
   pairs[2].column = eq_v1;

   for (int i = 0; i < 3; ++i)
     {
       matrix[row_non_zeros[eq_g] + i] = pairs[i].value;
       row_non_zeros[row_non_zeros[eq_g] + i] = pairs[i].column;
     }

   //Equation for V1
   solver->base_ws->matrix_rhs_storage->total_nz += 1;

   row_non_zeros[eq_v1 + 1] = row_non_zeros[eq_v1];

   rhs[eq_v1] = 0;

   matrix[eq_v1] = 1;

   //Equation for V2
   solver->base_ws->matrix_rhs_storage->total_nz += 1;

   row_non_zeros[eq_v2 + 1] = row_non_zeros[eq_v2];

   rhs[eq_v2] = 0;

   matrix[eq_v2] = 1;

   return 0;
}

int fill_dot_2 (solver_t *solver, int dot_number)
{
//  int top = solver->top_neighbor[dot_number];
  int bot = solver->bottom_neighbor[dot_number];
//  int left = solver->left_neighbor[dot_number];
//  int right = solver->right_neighbor[dot_number];

  int eq_g = get_equation_g (dot_number);
  int eq_v1 = get_equation_v1 (dot_number);
  int eq_v2 = get_equation_v2 (dot_number);

  double *prev = solver->base_ws->prev_answer;

  double *matrix = solver->base_ws->matrix_rhs_storage->matrix;
  double *rhs = solver->base_ws->matrix_rhs_storage->rhs;
  int *row_non_zeros = solver->base_ws->matrix_rhs_storage->row_non_zeros;
//  int n = solver->base_ws->matrix_rhs_storage->n;

  double tau = solver->tau;
//  double hx = solver->hx;
  double hy = solver->hy;

//  double mu = solver->mu;

//  double exp_minusG = exp (-prev[eq_g]);

  msr_pair pairs[8];

//  double mu_voln = prev[eq_g];

//  for (int i = 0; i < n; i += 3)
//    if (prev[i] < mu_voln)
//      mu_voln = prev[i];

//  mu_voln = mu * exp (-mu_voln);

  // Equation for G
   solver->base_ws->matrix_rhs_storage->total_nz += 4;

   rhs[eq_g] = 2 * prev[eq_g] + tau * prev[eq_g] * (prev[eq_v2] - prev[get_equation_v2 (bot)]) / hy + 2 * tau * f0 () -
       tau * (prev[eq_g] * prev[eq_v2] - 2.5 * prev[get_equation_g (bot)] * prev[get_equation_v2 (bot)]
       + 2 * prev[get_equation_g (solver->bottom_neighbor[bot])] * prev[get_equation_v2 (solver->bottom_neighbor[bot])] -
       0.5 * prev[get_equation_g (solver->bottom_neighbor[solver->bottom_neighbor[bot]])] * prev[get_equation_v2 (solver->bottom_neighbor[solver->bottom_neighbor[bot]])] +
       (2 - prev[eq_g]) * (prev[eq_v2] - 2.5 * prev[get_equation_v2 (bot)] + 2 * prev[get_equation_v2 (solver->bottom_neighbor[bot])] - 0.5 * prev[get_equation_v2 (solver->bottom_neighbor[solver->bottom_neighbor[bot]])])) / hy;

   row_non_zeros[eq_g + 1] = row_non_zeros[eq_g] + 3;

   matrix[eq_g] = 2 + tau * prev[eq_v2] / hy;

   pairs[0].value = -tau * prev[get_equation_v2 (bot)] / hy;
   pairs[0].column = get_equation_g (bot);

   pairs[1].value = -2 * tau / hy;
   pairs[1].column = get_equation_v2 (bot);

   pairs[2].value = 2 * tau / hy;
   pairs[2].column = eq_v2;

   for (int i = 0; i < 3; ++i)
     {
       matrix[row_non_zeros[eq_g] + i] = pairs[i].value;
       row_non_zeros[row_non_zeros[eq_g] + i] = pairs[i].column;
     }

   //Equation for V1
   solver->base_ws->matrix_rhs_storage->total_nz += 1;

   row_non_zeros[eq_v1 + 1] = row_non_zeros[eq_v1];

   rhs[eq_v1] = 0;

   matrix[eq_v1] = 1;

   //Equation for V2
   solver->base_ws->matrix_rhs_storage->total_nz += 1;

   row_non_zeros[eq_v2 + 1] = row_non_zeros[eq_v2];

   rhs[eq_v2] = 0;

   matrix[eq_v2] = 1;

   return 0;
}

int fill_dot_3 (solver_t *solver, int dot_number)
{
  int top = solver->top_neighbor[dot_number];
//  int bot = solver->bottom_neighbor[dot_number];
//  int left = solver->left_neighbor[dot_number];
//  int right = solver->right_neighbor[dot_number];

  int eq_g = get_equation_g (dot_number);
  int eq_v1 = get_equation_v1 (dot_number);
  int eq_v2 = get_equation_v2 (dot_number);

  double *prev = solver->base_ws->prev_answer;

  double *matrix = solver->base_ws->matrix_rhs_storage->matrix;
  double *rhs = solver->base_ws->matrix_rhs_storage->rhs;
  int *row_non_zeros = solver->base_ws->matrix_rhs_storage->row_non_zeros;
//  int n = solver->base_ws->matrix_rhs_storage->n;

  double tau = solver->tau;
//  double hx = solver->hx;
  double hy = solver->hy;

//  double mu = solver->mu;

//  double exp_minusG = exp (-prev[eq_g]);

  msr_pair pairs[8];

//  double mu_voln = prev[eq_g];

//  for (int i = 0; i < n; i += 3)
//    if (prev[i] < mu_voln)
//      mu_voln = prev[i];

//  mu_voln = mu * exp (-mu_voln);

  // Equation for G
   solver->base_ws->matrix_rhs_storage->total_nz += 4;

   rhs[eq_g] = 2 * prev[eq_g] + tau * (prev[get_equation_v2 (top)] - prev[eq_v2]) / hy + 2 * tau * f0 () +
       tau * (prev[eq_g] * prev[eq_v2] - 2.5 * prev[get_equation_g (top)] * prev[get_equation_v2 (top)]
       + 2 * prev[get_equation_g (solver->top_neighbor[top])] * prev[get_equation_v2 (solver->top_neighbor[top])] -
       0.5 * prev[get_equation_g (solver->top_neighbor[solver->top_neighbor[top]])] * prev[get_equation_v2 (solver->top_neighbor[solver->top_neighbor[top]])] +
       (2 - prev[eq_g]) * (prev[eq_v2] - 2.5 * prev[get_equation_v2 (top)] + 2 * prev[get_equation_v2 (solver->top_neighbor[top])] - 0.5 * prev[get_equation_v2 (solver->top_neighbor[solver->top_neighbor[top]])])) / hy;

   row_non_zeros[eq_g + 1] = row_non_zeros[eq_g] + 3;

   matrix[eq_g] = 2 - tau * prev[eq_v2] / hy;

   pairs[0].value = tau * prev[get_equation_v2 (top)] / hy;
   pairs[0].column = get_equation_g (top);

   pairs[1].value = 2 * tau / hy;
   pairs[1].column = get_equation_v2 (top);

   pairs[2].value = -2 * tau / hy;
   pairs[2].column = eq_v2;

   for (int i = 0; i < 3; ++i)
     {
       matrix[row_non_zeros[eq_g] + i] = pairs[i].value;
       row_non_zeros[row_non_zeros[eq_g] + i] = pairs[i].column;
     }

   //Equation for V1
   solver->base_ws->matrix_rhs_storage->total_nz += 1;

   row_non_zeros[eq_v1 + 1] = row_non_zeros[eq_v1];

   rhs[eq_v1] = 0;

   matrix[eq_v1] = 1;

   //Equation for V2
   solver->base_ws->matrix_rhs_storage->total_nz += 1;

   row_non_zeros[eq_v2 + 1] = row_non_zeros[eq_v2];

   rhs[eq_v2] = 0;

   matrix[eq_v2] = 1;

   return 0;
}

int fill_dot_6 (solver_t *solver, int dot_number)
{
  int top = solver->top_neighbor[dot_number];
//  int bot = solver->bottom_neighbor[dot_number];
//  int left = solver->left_neighbor[dot_number];
//  int right = solver->right_neighbor[dot_number];

  int eq_g = get_equation_g (dot_number);
  int eq_v1 = get_equation_v1 (dot_number);
  int eq_v2 = get_equation_v2 (dot_number);

  double *prev = solver->base_ws->prev_answer;

  double *matrix = solver->base_ws->matrix_rhs_storage->matrix;
  double *rhs = solver->base_ws->matrix_rhs_storage->rhs;
  int *row_non_zeros = solver->base_ws->matrix_rhs_storage->row_non_zeros;
//  int n = solver->base_ws->matrix_rhs_storage->n;

  double tau = solver->tau;
//  double hx = solver->hx;
  double hy = solver->hy;

//  double mu = solver->mu;

//  double exp_minusG = exp (-prev[eq_g]);

  msr_pair pairs[8];

//  double mu_voln = prev[eq_g];

//  for (int i = 0; i < n; i += 3)
//    if (prev[i] < mu_voln)
//      mu_voln = prev[i];

//  mu_voln = mu * exp (-mu_voln);

  // Equation for G
   solver->base_ws->matrix_rhs_storage->total_nz += 4;

   rhs[eq_g] = 2 * prev[eq_g] + tau * (prev[get_equation_v2 (top)] - prev[eq_v2]) / hy + 2 * tau * f0 () +
       tau * (prev[eq_g] * prev[eq_v2] - 2.5 * prev[get_equation_g (top)] * prev[get_equation_v2 (top)]
       + 2 * prev[get_equation_g (solver->top_neighbor[top])] * prev[get_equation_v2 (solver->top_neighbor[top])] -
       0.5 * prev[get_equation_g (solver->top_neighbor[solver->top_neighbor[top]])] * prev[get_equation_v2 (solver->top_neighbor[solver->top_neighbor[top]])] +
       (2 - prev[eq_g]) * (prev[eq_v2] - 2.5 * prev[get_equation_v2 (top)] + 2 * prev[get_equation_v2 (solver->top_neighbor[top])] - 0.5 * prev[get_equation_v2 (solver->top_neighbor[solver->top_neighbor[top]])])) / hy;

   row_non_zeros[eq_g + 1] = row_non_zeros[eq_g] + 3;

   matrix[eq_g] = 2 - tau * prev[eq_v2] / hy;

   pairs[0].value = tau * prev[get_equation_v2 (top)] / hy;
   pairs[0].column = get_equation_g (top);

   pairs[1].value = 2 * tau / hy;
   pairs[1].column = get_equation_v2 (top);

   pairs[2].value = -2 * tau / hy;
   pairs[2].column = eq_v2;

   for (int i = 0; i < 3; ++i)
     {
       matrix[row_non_zeros[eq_g] + i] = pairs[i].value;
       row_non_zeros[row_non_zeros[eq_g] + i] = pairs[i].column;
     }

   //Equation for V1
   solver->base_ws->matrix_rhs_storage->total_nz += 1;

   row_non_zeros[eq_v1 + 1] = row_non_zeros[eq_v1];

   rhs[eq_v1] = 0;

   matrix[eq_v1] = 1;

   //Equation for V2
   solver->base_ws->matrix_rhs_storage->total_nz += 2;

   row_non_zeros[eq_v2 + 1] = row_non_zeros[eq_v2] + 1;

   rhs[eq_v2] = 0;

   matrix[eq_v2] = -1;

   matrix[row_non_zeros[eq_v2]] = 1;
   row_non_zeros[row_non_zeros[eq_v2]] = get_equation_v2 (top);

   return 0;
}
