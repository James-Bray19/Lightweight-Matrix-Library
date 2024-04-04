/*
 * Lightweight Matrix Library (LML) - Source Code
 * File: lml.c
 * Author: James Bray
 *
 * This source file contains the implementation of the functions provided
 * by the Lightweight Matrix Library (LML). The library offers functionality
 * for generating matrices, retrieving data, performing matrix operations, editing
 * matrices, and various miscellaneous functions. The implementation aims for
 * efficiency and portability, making it suitable for use in embedded systems
 * and other resource-constrained environments.
 * 
 * For more detailed documentation, see the README in the following repo:
 *      https://github.com/James-Bray19/Lightweight-Matrix-Library
 * 
 * Any contributions or issues should also be posted there.
 */

#include "lml.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

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
        return 0.0;
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

Matrix *add(Matrix *mat1, Matrix *mat2) {
    // check if matrices have compatible dimensions
    if (mat1->rows != mat2->rows || mat1->cols != mat2->cols) {
        printf("Matrices aren't the same size.");
        return NULL;
    }

    // add matrices
    Matrix *result = zeros(mat1->rows, mat1->cols);
    for (int i = 0; i < result->rows; i++) {
        for (int j = 0; j < result->cols; j++) {
            result->data[i][j] = mat1->data[i][j] + mat2->data[i][j];
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

void *LU_decompose(Matrix *mat, Matrix **L, Matrix **U) {
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

Matrix *solve(Matrix *mat1, Matrix *mat2) {
    // check if the matrices are valid
    if (mat1 == NULL || mat2 == NULL) { printf("Invalid matrices\n"); return NULL; }

    // check if the coefficient matrix is square
    if (mat1->rows != mat1->cols) { printf("Coefficient matrix must be square\n"); return NULL; }

    // check if the dimensions match
    if (mat1->rows != mat2->rows) { printf("Dimension mismatch\n"); return NULL; }

    // perform LU decomposition
    Matrix *L, *U;
    LU_decompose(mat1, &L, &U);

    // solve Ly = mat2 using forward substitution
    Matrix *y = zeros(mat1->rows, 1);
    for (int i = 0; i < mat1->rows; i++) {
        double sum = 0.0;
        for (int j = 0; j < i; j++) {
            sum += L->data[i][j] * y->data[j][0];
        }
        y->data[i][0] = (mat2->data[i][0] - sum) / L->data[i][i];
    }

    // solve Ux = y using back substitution
    Matrix *x = zeros(mat1->rows, 1);
    for (int i = mat1->rows - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < mat1->rows; j++) {
            sum += U->data[i][j] * x->data[j][0];
        }
        x->data[i][0] = (y->data[i][0] - sum) / U->data[i][i];
    }

    release(L); release(U); release(y);

    return x;
}

Matrix *inverse(Matrix *mat) {
    // check if the matrix is square
    if (mat->rows != mat->cols) {
        printf("Matrix must be square\n");
        return NULL;
    }

    Matrix *I = identity(mat->rows);

    // perform LU decomposition
    Matrix *L = NULL;
    Matrix *U = NULL;
    LU_decompose(mat, &L, &U);

    Matrix *inverse_mat = zeros(mat->rows, mat->cols);

    // gaussian inverse method
    for (int i = 0; i < mat->cols; i++) {
        Matrix *column_i = get_col(I, i);
        Matrix *sol = solve(mat, column_i);
        set_col(inverse_mat, i, sol);

        release(sol); release(column_i);
    }

    release(L); release(U); release(I);
    return inverse_mat;
}

// --------------- In-Place Operations ---------------

void scale(Matrix *mat, double scalar) {
    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; j < mat->cols; j++) {
            mat->data[i][j] *= scalar;
        }   
    }
}

void shift(Matrix *mat, double scalar) {
    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; j < mat->cols; j++) {
            mat->data[i][j] += scalar;
        }   
    }
}

void map(Matrix *mat, double (*function)(double)) {
    // Check if the pointer to the matrix is valid
    if (mat == NULL) {
        printf("Invalid matrix pointer\n");
        return;
    }

    // Apply the function element-wise to the matrix
    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; j < mat->cols; j++) {
            mat->data[i][j] = function(mat->data[i][j]);
        }
    }
}

void set_row(Matrix *mat, int row_index, Matrix *row_values) {
    // check if row_index is valid
    if (row_index < 0 || row_index >= mat->rows) {
        printf("Invalid row index\n");
        return;
    }

    // check if dimensions are valid
    if (row_values->cols != mat->cols) {
        printf("Dimension mismatch\n");
        return;
    }

    // copy the values from row_values to the matrix
    for (int j = 0; j < mat->cols; j++) {
        mat->data[row_index][j] = row_values->data[0][j];
    }
}

void set_col(Matrix *mat, int col_index, Matrix *col_values) {
    // check if col_index is valid
    if (col_index < 0 || col_index >= mat->cols) {
        printf("Invalid column index\n");
        return;
    }

    // check if dimensions are valid
    if (col_values->rows != mat->rows) {
        printf("Dimension mismatch\n");
        return;
    }

    // copy the values from col_values to the matrix
    for (int i = 0; i < mat->rows; i++) {
        mat->data[i][col_index] = col_values->data[i][0];
    }
}

void remove_row(Matrix *mat, int row) {
    // check if row index is valid
    if (row < 0 || row >= mat->rows) {
        printf("Invalid row index\n");
        return;
    }

    // free memory for the row to be removed
    free(mat->data[row]);

    // shift rows above the removed row down
    for (int i = row + 1; i < mat->rows; i++) {
        mat->data[i - 1] = mat->data[i];
    }
    mat->rows--;

    // reallocate memory for the data array
    mat->data = realloc(mat->data, mat->rows * sizeof(double *));
}

void remove_col(Matrix *mat, int col) {
    // check if column index is valid
    if (col < 0 || col >= mat->cols) {
        printf("Invalid column index\n");
        return;
    }
    mat->cols--;

    // shift columns to the right of the removed column to the left
    for (int i = 0; i < mat->rows; i++) {
        for (int j = col; j < mat->cols; j++) {
            mat->data[i][j] = mat->data[i][j + 1];
        }
        // reallocate memory for the row to remove the column
        mat->data[i] = realloc(mat->data[i], mat->cols * sizeof(double));
    }
}

void insert_row(Matrix *mat, int row, Matrix *row_values) {
    // check pointer
    if (mat == NULL) { printf("Invalid matrix pointer\n"); return; }

    // check if row index is valid
    if (row < 0 || row > mat->rows) { printf("Invalid row index\n"); return; }

    // check if dimensions match
    if (row_values->cols != mat->cols) { printf("Dimension mismatch\n"); return; }

    mat->rows++;
    mat->data = realloc(mat->data, mat->rows * sizeof(double *));

    // shift existing rows down
    for (int i = mat->rows - 1; i > row; i--) {
        mat->data[i] = mat->data[i - 1];
    }

    // allocate memory for new row
    mat->data[row] = malloc(mat->cols * sizeof(double));

    // copy new row values
    for (int j = 0; j < mat->cols; j++) {
        mat->data[row][j] = row_values->data[0][j];
    }
}

void insert_col(Matrix *mat, int col, Matrix *col_values) {
    // check if the pointer to the matrix is valid
    if (mat == NULL) { printf("Invalid matrix pointer\n"); return; }

    // check if column index is valid
    if (col < 0 || col > mat->cols) { printf("Invalid column index\n"); return; }

    // check if dimensions match
    if (col_values->rows != mat->rows) { printf("Dimension mismatch\n"); return; }

    mat->cols++;
    double **new_data = malloc(mat->rows * sizeof(double *));
    
    // copy existing data and insert new column values
    for (int i = 0; i < mat->rows; i++) {
        new_data[i] = malloc(mat->cols * sizeof(double));

        // after insertion, k=j-1, so copying can continue
        for (int j = 0, k = 0; j < mat->cols; j++, k++) {
            if (j == col) {
                new_data[i][j] = col_values->data[i][0];
                k--;
            } else {
                new_data[i][j] = mat->data[i][k];
            }
        }
    }

    // free existing data
    for (int i = 0; i < mat->rows; i++) {
        free(mat->data[i]);
    }
    free(mat->data);

    mat->data = new_data;
}

void append_rows(Matrix *mat1, Matrix *mat2) {
    // check pointers
    if (mat1 == NULL || mat2 == NULL) { printf("Invalid matrix pointer\n"); return; }
    
    // check dimensions
    if (mat1->cols != mat2->cols) { printf("Number of columns does not match\n"); return; }
    
    // resize the matrix to accommodate additional rows
    mat1->data = realloc(mat1->data, (mat1->rows + mat2->rows) * sizeof(double *));
    
    // copy data from mat2 to mat1
    for (int i = 0; i < mat2->rows; i++) {
        mat1->data[mat1->rows + i] = malloc(mat1->cols * sizeof(double));
        for (int j = 0; j < mat1->cols; j++) {
            mat1->data[mat1->rows + i][j] = mat2->data[i][j];
        }
    }
    
    mat1->rows += mat2->rows;
    release(mat2);
}

void append_cols(Matrix *mat1, Matrix *mat2) {
    // check pointers
    if (mat1 == NULL || mat2 == NULL) { printf("Invalid matrix pointer\n"); return; }
    
    // check dimensions
    if (mat1->rows != mat2->rows) { printf("Number of rows does not match\n"); return; }
    
    // resize the matrix to accommodate additional columns
    for (int i = 0; i < mat1->rows; i++) {
        mat1->data[i] = realloc(mat1->data[i], (mat1->cols + mat2->cols) * sizeof(double));
        for (int j = 0; j < mat2->cols; j++) {
            mat1->data[i][mat1->cols + j] = mat2->data[i][j];
        }
    }
    
    mat1->cols += mat2->cols;
    release(mat2);
}

// --------------- Miscellaneous Functions ---------------

void display(Matrix *mat) {
    if (mat == NULL) {
        printf("Matrix is NULL\n");
        return;
    }

    printf("Matrix (%d x %d):\n", mat->rows, mat->cols);

    // print matrix
    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; j < mat->cols; j++) {
            printf("%8.2f ", mat->data[i][j]);
        }
        printf("\n");
    }
}

void release(Matrix *mat) {
    if (mat != NULL) {
        // free memory for each row
        for (int i = 0; i < mat->rows; i++) {
            free(mat->data[i]);
        }
        
        free(mat->data);
        free(mat);
    }
}