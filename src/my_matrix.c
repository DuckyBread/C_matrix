#include "my_matrix.h"

/** Main Functions */
int my_create_matrix(int rows, int columns, matrix_t *result) {
  result->matrix = NULL;
  result->rows = 0;
  result->columns = 0;
  int res = OK;
  if (rows > 0 && columns > 0) {
    result->matrix = (double **)calloc((size_t)rows, (size_t)sizeof(double *));
    for (int c = 0; c < rows; c++) {
      result->matrix[c] =
          (double *)calloc((size_t)columns, (size_t)sizeof(double));
    }
    result->rows = rows;
    result->columns = columns;
  } else {
    res = ERROR_MATRIX;
  }
  return res;
}

void my_remove_matrix(matrix_t *A) {
  if (A->matrix) {
    for (int r = 0; r < A->rows; r++) {
      free(A->matrix[r]);
    }
    free(A->matrix);
  }
  A->columns = 0;
  A->rows = 0;
  A->matrix = NULL;
}

int my_eq_matrix(matrix_t *A, matrix_t *B) {
  int res = SUCCESS;

  if (A->rows == B->rows && A->columns == B->columns && A->matrix != NULL &&
      B->matrix != NULL) {
    for (int r = 0; r < A->rows; r++) {
      for (int c = 0; c < A->columns; c++) {
        if (fabs(*(*(A->matrix + r) + c) - *(*(B->matrix + r) + c)) >= 1E-6) {
          res = FAILURE;
          break;
        }
      }
    }
  } else {
    res = FAILURE;
  }
  return res;
}

int my_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res = OK;
  if (B->columns > 0 && B->rows > 0 && A->columns > 0 && A->rows > 0 &&
      A->matrix != NULL && B->matrix != NULL) {
    if (A->rows == B->rows && A->columns == B->columns) {
      my_create_matrix(A->rows, A->columns, result);
      for (int r = 0; r < A->rows; r++) {
        for (int c = 0; c < A->columns; c++) {
          *(*(result->matrix + r) + c) =
              *(*(A->matrix + r) + c) + *(*(B->matrix + r) + c);
        }
      }
    } else {
      res = ERROR_CALCULATION;
    }
  } else {
    res = ERROR_MATRIX;
  }

  return res;
}

int my_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res = OK;
  result->matrix = NULL;
  result->columns = 0;
  result->rows = 0;
  if (A->matrix == NULL || B->matrix == NULL) {
    res = ERROR_MATRIX;
  } else if (A->rows == B->rows && A->columns == B->columns && A->columns > 0 &&
             A->rows > 0) {
    my_create_matrix(A->rows, A->columns, result);
    for (int r = 0; r < A->rows; r++) {
      for (int c = 0; c < A->columns; c++) {
        *(*(result->matrix + r) + c) =
            *(*(A->matrix + r) + c) - *(*(B->matrix + r) + c);
      }
    }
  } else {
    res = ERROR_CALCULATION;
  }
  return res;
}

int my_mult_number(matrix_t *A, double number, matrix_t *result) {
  int res = OK;
  result->matrix = NULL;
  result->columns = 0;
  result->rows = 0;
  if (A->matrix == NULL) {
    res = ERROR_MATRIX;
  } else if (A->columns > 0 && A->rows > 0) {
    my_create_matrix(A->rows, A->columns, result);
    for (int r = 0; r < A->rows; r++) {
      for (int c = 0; c < A->columns; c++) {
        *(*(result->matrix + r) + c) = *(*(A->matrix + r) + c) * number;
      }
    }
  } else {
    res = ERROR_CALCULATION;
  }
  return res;
}

int my_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res = OK;
  result->matrix = NULL;
  result->columns = 0;
  result->rows = 0;
  if (A->matrix == NULL || B->matrix == NULL) {
    res = ERROR_MATRIX;
  } else if (A->columns == B->rows && A->columns > 0 && A->rows > 0) {
    my_create_matrix(A->rows, B->columns, result);

    for (int r = 0; r < result->rows; r++) {
      for (int c = 0; c < result->columns; c++) {
        double sum = 0;
        for (int c_A = 0, r_B = 0; c_A < A->columns; c_A++, r_B++) {
          sum += *(*(A->matrix + r) + c_A) * *(*(B->matrix + r_B) + c);
        }
        *(*(result->matrix + r) + c) = sum;
      }
    }

  } else {
    res = ERROR_CALCULATION;
  }
  return res;
}

int my_transpose(matrix_t *A, matrix_t *result) {
  int res = OK;
  result->matrix = NULL;
  result->columns = 0;
  result->rows = 0;
  if (A->matrix == NULL) {
    res = ERROR_MATRIX;
  } else if (A->columns > 0 && A->rows > 0) {
    my_create_matrix(A->columns, A->rows, result);

    for (int r = 0; r < result->rows; r++) {
      for (int c = 0; c < result->columns; c++) {
        *(*(result->matrix + r) + c) = *(*(A->matrix + c) + r);
      }
    }

  } else {
    res = ERROR_CALCULATION;
  }
  return res;
}

int my_calc_complements(matrix_t *A, matrix_t *result) {
  int res = OK;
  if (A->rows > 0 && A->columns > 0 && A->matrix != NULL) {
    int m_res = my_create_matrix(A->rows, A->columns, result);
    if (A->rows == A->columns && m_res == OK) {
      if (A->rows == 1) {
        *(*(result->matrix)) = 1;
      } else {
        for (int r = 0; r < A->rows; r++) {
          for (int c = 0; c < A->columns; c++) {
            matrix_t tmp;
            my_create_matrix(A->rows - 1, A->columns - 1, &tmp);
            minor_new(A, r, c, &tmp);
            *(*(result->matrix + r) + c) = determ(&tmp) * pow(-1.0, r + c);
            my_remove_matrix(&tmp);
          }
        }
      }
    } else {
      res = ERROR_CALCULATION;
    }
  } else {
    res = ERROR_MATRIX;
  }

  return res;
}

int my_determinant(matrix_t *A, double *result) {
  int res = OK;
  *result = NAN;
  if (A->rows > 0 && A->columns > 0 && A->matrix != NULL) {
    if (A->rows == A->columns) {
      *result = determ(A);
    } else {
      res = ERROR_CALCULATION;
    }
  } else {
    res = ERROR_MATRIX;
  }
  return res;
}

int my_inverse_matrix(matrix_t *A, matrix_t *result) {
  int res = OK;
  if (A->rows > 0 && A->columns > 0 && A->matrix != NULL) {
    double det_A = 0;
    int check = my_determinant(A, &det_A);
    if (A->rows == A->columns && check == OK && fabs(det_A) > 1e-6) {
      my_create_matrix(A->rows, A->columns, result);

      matrix_t alg_dop, transp;
      my_calc_complements(A, &alg_dop);
      my_transpose(&alg_dop, &transp);

      for (int r = 0; r < transp.rows; r++) {
        for (int c = 0; c < transp.rows; c++) {
          *(*(result->matrix + r) + c) = *(*(transp.matrix + r) + c) / det_A;
        }
      }

      my_remove_matrix(&alg_dop);
      my_remove_matrix(&transp);
    } else {
      res = ERROR_CALCULATION;
    }
  } else {
    res = ERROR_MATRIX;
  }
  return res;
}

/** Help Functions */

void my_print_matrix(matrix_t A) {
  for (int r = 0; r < A.rows; r++) {
    for (int c = 0; c < A.columns; c++) {
      printf("%f\t", *(*(A.matrix + r) + c));
    }
    printf("\n");
  }
}

void minor_new(matrix_t *A, int R, int C, matrix_t *tmp) {
  int tmp_r = 0;
  for (int r = 0; r < A->rows; r++) {
    int tmp_c = 0;
    if (R != r) {
      for (int c = 0; c < A->columns; c++) {
        if (C != c) {
          *(*(tmp->matrix + tmp_r) + tmp_c) = *(*(A->matrix + r) + c);
          tmp_c++;
        }
      }
      tmp_r++;
    }
  }
}

double determ(matrix_t *A) {
  double res = 0;
  if (A->rows == 1) {
    res = *(*(A->matrix));
  } else if (A->rows == 2) {
    res = *(*(A->matrix)) * *(*(A->matrix + 1) + 1) -
          *(*(A->matrix) + 1) * *(*(A->matrix + 1));
  } else {
    double sign = 1.;
    for (int col = 0; col < A->columns; col++) {
      matrix_t tmp;
      my_create_matrix(A->rows - 1, A->columns - 1, &tmp);
      minor_new(A, 0, col, &tmp);
      res += sign * *(*(A->matrix) + col) * determ(&tmp);
      sign *= -1;
      my_remove_matrix(&tmp);
    }
  }
  return res;
}
