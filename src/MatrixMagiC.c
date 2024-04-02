// matrixmagic.c

#include "matrixmagic.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

// --------------- Generating Matrices ---------------

Matrix *zeros(int rows, int cols) {

    // allocate memory for array
    Matrix *mat = malloc(sizeof(Matrix));
    if (mat == NULL) { return NULL; }
    mat->data = malloc(rows * sizeof(double *));
    if (mat->data == NULL) { free(mat); return NULL; }

    // initialize matrix elements to zero
    for (int i = 0; i < rows; i++) {
        mat->data[i] = calloc(cols, sizeof(double));

        // handle memory allocation failure
        if (mat->data[i] == NULL) {
            for (int j = 0; j < i; j++) { free(mat->data[j]); }
            free(mat->data); free(mat); return NULL;
        }
    }

    mat->rows = rows;
    mat->cols = cols;

    return mat;
}

Matrix *ones(int rows, int cols) {
    // allocate memory for array
    Matrix *mat = malloc(sizeof(Matrix));
    if (mat == NULL) { return NULL; }
    mat->data = malloc(rows * sizeof(double *));
    if (mat->data == NULL) { free(mat); return NULL; }

    // initialize matrix elements to one
    for (int i = 0; i < rows; i++) {
        mat->data[i] = calloc(cols, sizeof(double));

        // handle memory allocation failure
        if (mat->data[i] == NULL) {
            for (int j = 0; j < i; j++) { free(mat->data[j]); }
            free(mat->data); free(mat); return NULL;
        }

        for (int j = 0; j < cols; j++) { mat->data[i][j] = 1.0; }
    }
    mat->rows = rows;
    mat->cols = cols;

    return mat;
}

Matrix *identity(int size) {
    // allocate memory for the identity matrix
    Matrix *mat = zeros(size, size);

    // set leading diagonal to 1
    for (int i = 0; i < size; i++) { mat->data[i][i] = 1.0; }

    mat->rows = size;
    mat->cols = size;

    return mat;
}

Matrix *random(int rows, int cols) {

    // set generator seed
    srand((unsigned int)(time(NULL) * clock()));

    // allocate memory for array
    Matrix *mat = malloc(sizeof(Matrix));
    if (mat == NULL) { return NULL; }
    mat->data = malloc(rows * sizeof(double *));
    if (mat->data == NULL) { free(mat); return NULL; }

    // initialize matrix elements to random numbers 0-1
    for (int i = 0; i < rows; i++) {
        mat->data[i] = calloc(cols, sizeof(double));

        // handle memory allocation failure
        if (mat->data[i] == NULL) {
            for (int j = 0; j < i; j++) { free(mat->data[j]); }
            free(mat->data); free(mat); return NULL;
        }

        for (int j = 0; j < cols; j++) { mat->data[i][j] = (double)rand() / RAND_MAX; }
    }
    mat->rows = rows;
    mat->cols = cols;

    return mat;
}

Matrix *matrix_from_array(int rows, int cols, double array[rows][cols]) {
    // allocate memory for array
    Matrix *mat = malloc(sizeof(Matrix));
    if (mat == NULL) { return NULL; }
    mat->data = malloc(rows * sizeof(double *));
    if (mat->data == NULL) { free(mat); return NULL; }

    // initialize matrix elements to array elements
    for (int i = 0; i < rows; i++) {
        mat->data[i] = calloc(cols, sizeof(double));

        // handle memory allocation failure
        if (mat->data[i] == NULL) {
            for (int j = 0; j < i; j++) { free(mat->data[j]); }
            free(mat->data); free(mat); return NULL;
        }

        for (int j = 0; j < cols; j++) { mat->data[i][j] = array[i][j]; }
    }

    mat->rows = rows;
    mat->cols = cols;

    return mat;
}

// --------------- Retrieving Data ---------------

Matrix *get_row(Matrix *mat, int row) {
    // check if the row index is valid
    if (row < 0 || row >= mat->rows) {
        printf("Invalid row index\n");
        return NULL;
    }

    // create a new matrix to represent the row
    Matrix *row_matrix = zeros(1, mat->cols);

    // copy the row data into the new matrix
    for (int j = 0; j < mat->cols; j++) {
        row_matrix->data[0][j] = mat->data[row][j];
    }

    return row_matrix;
}

Matrix *get_col(Matrix *mat, int col) {
    // check if the column index is valid
    if (col < 0 || col >= mat->cols) {
        fprintf(stderr, "Invalid column index\n");
        return NULL;
    }

    // create a new matrix to represent the column
    Matrix *col_matrix = zeros(mat->rows, 1);

    // copy the column data into the new matrix
    for (int i = 0; i < mat->rows; i++) {
        col_matrix->data[i][0] = mat->data[i][col];
    }

    return col_matrix;
}

Matrix *copy(Matrix *mat) {

    // create a new matrix with the same dimensions as the original
    Matrix *copied_matrix = zeros(mat->rows, mat->cols);

    // copy the data from the original matrix to the new matrix
    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; j < mat->cols; j++) {
            copied_matrix->data[i][j] = mat->data[i][j];
        }
    }

    return copied_matrix;
}

Matrix *get_lower(Matrix *mat) {
    // create a new matrix
    Matrix *lower_triangular = zeros(mat->rows, mat->cols);

    // copy the lower triangular part
    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; j < mat->cols; j++) {
            if (i >= j) { lower_triangular->data[i][j] = mat->data[i][j]; }
        }
    }

    return lower_triangular;
}

Matrix *get_upper(Matrix *mat) {
    // create a new matrix
    Matrix *upper_triangular = zeros(mat->rows, mat->cols);

    // copy the upper triangular part
    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; j < mat->cols; j++) {
            if (i <= j) { upper_triangular->data[i][j] = mat->data[i][j]; }
        }
    }

    return upper_triangular;
}

Matrix *get_submatrix(Matrix *mat, int row, int col, int rows, int cols) {
    // check if the indices are within bounds
    if (row < 0 || col < 0 || row + rows > mat->rows || col + cols > mat->cols) {
        fprintf(stderr, "Submatrix indices out of bounds\n");
        return NULL;
    }

    // create a new matrix for the submatrix
    Matrix *submatrix = zeros(rows, cols);
    if (submatrix == NULL) {
        fprintf(stderr, "Failed to create submatrix\n");
        return NULL;
    }

    // copy the elements from the original matrix to the submatrix
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            submatrix->data[i][j] = mat->data[row + i][col + j];
        }
    }

    return submatrix;
}

// --------------- Matrix Operations ---------------

double det(Matrix *mat) {
    if (mat->rows != mat->cols) {
        printf("Determinant can only be calculated for square matrices\n");
        return NAN;
    }

    // perform LU decomposition
    Matrix *L, *U;
    LU_decompose(mat, &L, &U);

    // determinant is product of upper diagonal
    double det = 1.0;
    for (int i = 0; i < mat->rows; i++) { det *= U->data[i][i]; }

    release(L); release(U);
    return det;
}

Matrix *transpose(Matrix *mat) {
    // create a new matrix with dimensions swapped
    Matrix *transposed = zeros(mat->cols, mat->rows);
    if (transposed == NULL) {
        printf("Failed to create transposed matrix\n");
        return NULL;
    }

    // copy elements from the original matrix to the transposed matrix
    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; j < mat->cols; j++) {
            transposed->data[j][i] = mat->data[i][j];
        }
    }

    return transposed;
}

Matrix *matrix_inverse(Matrix *mat);

Matrix *scale(Matrix *mat, double scalar) {
    Matrix *result = zeros(mat->rows, mat->cols);

    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; j < mat->cols; j++) {
            result->data[i][j] = scalar * mat->data[i][j];
        }   
    }

    return result;
}

Matrix *multiply(Matrix *mat1, Matrix *mat2) {
    // check dimensions
    if (mat1->cols != mat2->rows) {
        printf("Matrices are not compatible for multiplication\n");
        return NULL;
    }

    // multiply matrices together
    Matrix *result = zeros(mat1->rows, mat2->cols);
    for (int i = 0; i < mat1->rows; i++) {
        for (int j = 0; j < mat2->cols; j++) {
            for (int k = 0; k < mat1->cols; k++) {
                result->data[i][j] += mat1->data[i][k] * mat2->data[k][j];
            }
        }
    }

    return result;
}

Matrix *LU_decompose(Matrix *mat, Matrix **L, Matrix **U) {
    if (mat->rows != mat->cols) {
        printf("LU decomposition requires a square matrix\n");
    }

    int n = mat->rows;

    // initialise
    *L = identity(n);
    *U = zeros(n, n);

    // perform LU decomposition
    for (int j = 0; j < n; j++) {
        (*U)->data[0][j] = mat->data[0][j];
    }
    for (int i = 1; i < n; i++) {
        (*L)->data[i][0] = mat->data[i][0] / (*U)->data[0][0];
    }
    for (int i = 1; i < n; i++) {
        for (int j = i; j < n; j++) {
            double sum = 0.0;
            for (int k = 0; k < i; k++) {
                sum += (*L)->data[i][k] * (*U)->data[k][j];
            }
            (*U)->data[i][j] = mat->data[i][j] - sum;
        }
        for (int j = i + 1; j < n; j++) {
            double sum = 0.0;
            for (int k = 0; k < i; k++) {
                sum += (*L)->data[j][k] * (*U)->data[k][i];
            }
            (*L)->data[j][i] = (mat->data[j][i] - sum) / (*U)->data[i][i];
        }
    }
}

Matrix *solve_system(Matrix *mat1, Matrix *mat2);


// --------------- Matrix Editing ---------------

Matrix *map(Matrix *mat, double (*function)(double));
void set_row(Matrix *mat, int row_index, double *row_values);
void set_col(Matrix *mat, int col_index, double *col_values);
Matrix *remove_rows(Matrix *mat, int start_row, int num_rows);
Matrix *remove_cols(Matrix *mat, int start_col, int num_cols);
Matrix *insert_row(Matrix *mat, int row_index, double *row_values);
Matrix *insert_col(Matrix *mat, int col_index, double *col_values);

// --------------- Miscellaneous Functions ---------------

void display(Matrix *mat) {
    if (mat == NULL) {
        printf("Matrix is NULL\n");
        return;
    }

    printf("Matrix (%d x %d):\n", mat->rows, mat->cols);
    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; j < mat->cols; j++) {
            printf("%8.2f ", mat->data[i][j]);
        }
        printf("\n");
    }
}

void release(Matrix *mat) {
    if (mat != NULL) {
        // Free memory for each row
        for (int i = 0; i < mat->rows; i++) {
            free(mat->data[i]);
        }
        // Free memory for the array of row pointers
        free(mat->data);
        // Free memory for the Matrix structure itself
        free(mat);
    }
}