#ifndef _SRC_MY_MATRIX_H_
#define _SRC_MY_MATRIX_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/** For Equalization */
#define SUCCESS 1
#define FAILURE 0
/** ---------------- */

#define OK 0
/** INCORRECT MATRIX */
#define ERROR_MATRIX 1
/** MISMATCHED MATRIX SIZES
 *  MATRIX FOR WHICH CALCULATIONS CAN NOT BE PERFORMED */
#define ERROR_CALCULATION 2

typedef struct matrix_struct {
  double **matrix;
  int rows;
  int columns;
} matrix_t;

/** Main Functions */

int my_create_matrix(int rows, int columns, matrix_t *result);
void my_remove_matrix(matrix_t *A);
int my_eq_matrix(matrix_t *A, matrix_t *B);
int my_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int my_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int my_mult_number(matrix_t *A, double number, matrix_t *result);
int my_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int my_transpose(matrix_t *A, matrix_t *result);
int my_calc_complements(matrix_t *A, matrix_t *result);
int my_determinant(matrix_t *A, double *result);
int my_inverse_matrix(matrix_t *A, matrix_t *result);

/** Help Functions */

void my_print_matrix(matrix_t A);
// void my_minor(matrix_t origin, int orig_r, int orig_col, matrix_t *tmp);
// int my_det(matrix_t A, double *result);
// int my_minor_matrix(matrix_t origin, matrix_t *result);

void minor_new(matrix_t *A, int R, int C, matrix_t *tmp);
double determ(matrix_t *A);

#endif  // SRC_MY_MATRIX_H_