#ifndef SOLVER
#define SOLVER

typedef struct
{
  int n;
  int total_nz;
  double *matrix;
  double *rhs;
  int *row_non_zeros;
} msr_matrix_rhs;

typedef struct
{
  int iter;
  int max_iter;

  int total_dots_amount;
  int total_equation_amount;


  msr_matrix_rhs *matrix_rhs_storage;
  msr_matrix_rhs *matrix_rhs_storage_preconditioned;

  double *prev_answer;

  double *r_start;
  double *r;
  double *p;
  double *u;
  double *q;
  double *tmp1, tmp2;

  double alpha;
  double beta;
  double last_residual;
} solver_ws;

typedef struct
{
  int MX;
  int MY;
  int N;
  double T;
  double hx;
  double hy;
  double tau;

  double border_omega;
  double mu;

  int *top_neighbor;
  int *bottom_neighbor;
  int *left_neighbor;
  int *right_neighbor;
  int *dot_type;

  /*
   *       ___2___
   *   1->|       |
   *      |___  0 |
   *        3 |   |<-5
   *          |   |
   *       4->|   |
   *          |___|
   *            6
   */

  // Если не будет сходиться, то можно исправить угловые точки

  int solver_type; // 1 - Laspack, 2 - own BiCGS

  solver_ws *base_ws;
} solver_t;

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
                 int solver_type_arg);

void solver_free (solver_t *solver);

void solver_run (solver_t *solver);

void msr_matrix_rhs_init (int total_equations_amount,
                          msr_matrix_rhs *matrix_rhs_storage);

int get_equation_g (int dot_number);
int get_equation_v1 (int dot_number);
int get_equation_v2 (int dot_number);




#endif
