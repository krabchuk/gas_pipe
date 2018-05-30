#include "neighbors_and_type_init.h"
#include <stdio.h>
#include <stdlib.h>

void init_neighbors_and_type (solver_t *solver)
{
  int n = solver->base_ws->total_dots_amount;

  int MX = solver->MX;
  int MY = solver->MY;

  int last_row;
  int first_num;
  int first_num_in_row;
  int first_num_in_last_row;

  solver->top_neighbor = (int *)malloc(n * sizeof (int));
  solver->bottom_neighbor = (int *)malloc(n * sizeof (int));
  solver->left_neighbor = (int *)malloc(n * sizeof (int));
  solver->right_neighbor = (int *)malloc(n * sizeof (int));
  solver->dot_type = (int *)malloc(n * sizeof (int));

  // Bottom border

  solver->top_neighbor[0] = MX + 1;
  solver->bottom_neighbor[0] = -1;
  solver->left_neighbor[0] = -1;
  solver->right_neighbor[0] = 1;
  solver->dot_type[0] = 4;

  for (int i = 1; i < MX; ++i)
    {
      solver->top_neighbor[i] = i + MX + 1;
      solver->bottom_neighbor[i] = -1;
      solver->left_neighbor[i] = i - 1;
      solver->right_neighbor[i] = i + 1;
      solver->dot_type[i] = 6;
    }

  solver->top_neighbor[MX] = MX + MX +1;
  solver->bottom_neighbor[MX] = -1;
  solver->left_neighbor[MX] = MX - 1;
  solver->right_neighbor[MX] = -1;
  solver->dot_type[MX] = 5;


  for (int row = 1; row < MY * 2 - 1; ++row)
    {
      int first_num_in_row = row * MX + 1;

      // Left border
      solver->top_neighbor[first_num_in_row] = first_num_in_row + MX + 1;
      solver->bottom_neighbor[first_num_in_row] = first_num_in_row - MX - 1;
      solver->left_neighbor[0] = -1;
      solver->right_neighbor[0] = first_num_in_row + 1;
      solver->dot_type[0] = 4;

      // Internal dots
      for (int i = first_num_in_row; i < first_num_in_row + MX; ++i)
        {
          solver->top_neighbor[i] = i + MX + 1;
          solver->bottom_neighbor[i] = i - MX - 1;
          solver->left_neighbor[i] = i - 1;
          solver->right_neighbor[i] = i + 1;
          solver->dot_type[i] = 0;
        }

      // Right border
      solver->top_neighbor[first_num_in_row + MX] = first_num_in_row + MX + MX + 1;
      solver->bottom_neighbor[first_num_in_row + MX] = first_num_in_row - 1;
      solver->left_neighbor[first_num_in_row + MX] = first_num_in_row + MX - 1;
      solver->right_neighbor[first_num_in_row + MX] = -1;
      solver->dot_type[MX] = 5;
    }

  last_row = MY * 2 - 1;
  first_num_in_row = last_row * MX + 1;

  // Left border
  solver->top_neighbor[first_num_in_row] = first_num_in_row + 2 * MX + 1;
  solver->bottom_neighbor[first_num_in_row] = first_num_in_row - MX - 1;
  solver->left_neighbor[0] = -1;
  solver->right_neighbor[0] = first_num_in_row + 1;
  solver->dot_type[0] = 4;

  // Internal dots
  for (int i = first_num_in_row; i < first_num_in_row + MX; ++i)
    {
      solver->top_neighbor[i] = i + 2 * MX + 1;
      solver->bottom_neighbor[i] = i - MX - 1;
      solver->left_neighbor[i] = i - 1;
      solver->right_neighbor[i] = i + 1;
      solver->dot_type[i] = 0;
    }

  // Right border
  solver->top_neighbor[first_num_in_row + MX] = first_num_in_row + MX + 2 * MX + 1;
  solver->bottom_neighbor[first_num_in_row + MX] = first_num_in_row - 1;
  solver->left_neighbor[first_num_in_row + MX] = first_num_in_row + MX - 1;
  solver->right_neighbor[first_num_in_row + MX] = -1;
  solver->dot_type[MX] = 5;


  /*-------------------End of vertical part--------------------*/
  /*-------------------Start of horizontal part----------------*/

  first_num = 2 * MY * (MX + 1);

  solver->top_neighbor[first_num] = first_num + 2 * MX + 1;
  solver->bottom_neighbor[first_num] = -1;
  solver->left_neighbor[first_num] = -1;
  solver->right_neighbor[first_num] = first_num + 1;
  solver->dot_type[first_num] = 3;

  for (int i = first_num + 1; i < first_num + MX; ++i)
    {
      solver->top_neighbor[i] = i + 2 * MX + 1;
      solver->bottom_neighbor[i] = -1;
      solver->left_neighbor[i] = i - 1;
      solver->right_neighbor[i] = i + 1;
      solver->dot_type[i] = 3;
    }

  solver->top_neighbor[first_num + MX] = first_num + MX + 2 * MX + 1;
  solver->bottom_neighbor[first_num + MX] = first_num + MX - 2 * MX - 1;
  solver->left_neighbor[first_num + MX] = first_num + MX - 1;
  solver->right_neighbor[first_num + MX] = first_num + MX + 1;
  solver->dot_type[first_num + MX] = 3;

  for (int i = first_num + MX + 1; i < first_num + 2 * MX; ++i)
    {
      solver->top_neighbor[i] = i + 2 * MX + 1;
      solver->bottom_neighbor[i] = i - 2 * MX - 1;
      solver->left_neighbor[i] = i - 1;
      solver->right_neighbor[i] = i + 1;
      solver->dot_type[i] = 0;
    }

  solver->top_neighbor[first_num + 2 * MX] = first_num + 2 * MX + 2 * MX + 1;
  solver->bottom_neighbor[first_num + 2 * MX] = first_num + 2 * MX - 2 * MX - 1;
  solver->left_neighbor[first_num + 2 * MX] = first_num + 2 * MX - 1;
  solver->right_neighbor[first_num + 2 * MX] = -1;
  solver->dot_type[first_num + 2 * MX] = 5;

  for (int row = 1; row < MY; ++row)
    {
      first_num_in_row = first_num + 2 * MX * row + 1;

      solver->top_neighbor[first_num_in_row] = first_num_in_row + 2 * MX + 1;
      solver->bottom_neighbor[first_num_in_row] = first_num_in_row - 2 * MX - 1;
      solver->left_neighbor[first_num_in_row] = -1;
      solver->right_neighbor[first_num_in_row] = first_num_in_row + 1;
      solver->dot_type[first_num_in_row] = 1;

      for (int i = first_num_in_row; i < first_num_in_row + 2 * MX; ++i)
        {
          solver->top_neighbor[i] = i + 2 * MX + 1;
          solver->bottom_neighbor[i] = i - 2 * MX - 1;
          solver->left_neighbor[i] = i - 1;
          solver->right_neighbor[i] = i + 1;
          solver->dot_type[i] = 0;
        }

      solver->top_neighbor[first_num_in_row + 2 * MX] = first_num_in_row + 2 * MX + 2 * MX + 1;
      solver->bottom_neighbor[first_num_in_row + 2 * MX] = first_num_in_row + 2 * MX - 2 * MX - 1;
      solver->left_neighbor[first_num_in_row + 2 * MX] = first_num_in_row + 2 * MX - 1;
      solver->right_neighbor[first_num_in_row + 2 * MX] = -1;
      solver->dot_type[first_num_in_row + 2 * MX] = 5;
    }

  first_num_in_last_row = first_num + (2 * MX + 1) * MY;

  solver->top_neighbor[first_num_in_last_row] = -1;
  solver->bottom_neighbor[first_num_in_last_row] = first_num_in_last_row - 2 * MX - 1;
  solver->left_neighbor[first_num_in_last_row] = -1;
  solver->right_neighbor[first_num_in_last_row] = first_num_in_last_row + 1;
  solver->dot_type[first_num_in_last_row] = 2;

  for (int i = first_num_in_last_row; i < first_num_in_last_row + 2 * MX; ++i)
    {
      solver->top_neighbor[i] = -1;
      solver->bottom_neighbor[i] = i - 2 * MX - 1;
      solver->left_neighbor[i] = i - 1;
      solver->right_neighbor[i] = i + 1;
      solver->dot_type[i] = 2;
    }

  solver->top_neighbor[first_num_in_last_row + 2 * MX] = -1;
  solver->bottom_neighbor[first_num_in_last_row + 2 * MX] = first_num_in_last_row + 2 * MX - 2 * MX - 1;
  solver->left_neighbor[first_num_in_last_row + 2 * MX] = first_num_in_last_row + 2 * MX - 1;
  solver->right_neighbor[first_num_in_last_row + 2 * MX] = -1;
  solver->dot_type[first_num_in_last_row + 2 * MX] = 2;

}
