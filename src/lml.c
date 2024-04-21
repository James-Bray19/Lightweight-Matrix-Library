/*
 * Lightweight Matrix Library (LML) - Source Code
 *
 * File:     lml.c
 * Author:   James Bray
 * Repo:     https://github.com/James-Bray19/Lightweight-Matrix-Library
 * 
 * Implementation file for the LML (Matrix Library) containing
 * definitions of functions declared in the lml.h header file.
 */

#include "lml.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

// --------------- Static Functions ---------------

static void lml_error(const char *function_name, const char *error) {
    // format the error
    printf("\n[LML Error]: in function '%s' - %s.\n", function_name, error);
}

// --------------- Generating Matrices ---------------

Matrix *zeros(int rows, int cols) {
    
    // input validation
    if (rows <= 0 || cols <= 0) { lml_error(__func__, "invalid dimensions"); return NULL; }

    // allocate memory for matrix
    Matrix *mat = malloc(sizeof(Matrix));
    if (mat == NULL) {
        lml_error(__func__, "memory allocation failed");
        return NULL; 
    }

    // allocate memory for data
    mat->data = malloc(rows * sizeof(double *));
    if (mat->data == NULL) { 
        free(mat); 
        lml_error(__func__, "memory allocation failed");
        return NULL; 
    }

    // initialize matrix elements to zero
    for (int i = 0; i < rows; i++) {
        mat->data[i] = calloc(cols, sizeof(double));

        // handle errors
        if (mat->data[i] == NULL) {
            for (int j = 0; j < i; j++) { free(mat->data[j]); }
            free(mat->data); free(mat); 
            lml_error(__func__, "memory allocation failed");
            return NULL;
        }
    }

    mat->rows = rows;
    mat->cols = cols;

    return mat;
}

Matrix *constants(int rows, int cols, int value) {
    
    // input validation
    if (rows <= 0 || cols <= 0) { lml_error(__func__, "invalid dimensions"); return NULL; }

    // allocate memory for matrix
    Matrix *mat = zeros(rows, cols);

    // initialize matrix elements to specific value
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) { 
            mat->data[i][j] = value; 
        }
    }

    return mat;
}

Matrix *identity(int size) {
    
    // input validation
    if (size <= 0) { lml_error(__func__, "invalid size"); return NULL; }

    // allocate memory for matrix
    Matrix *mat = zeros(size, size);

    // set leading diagonal to 1
    for (int i = 0; i < size; i++) { mat->data[i][i] = 1.0; }

    mat->rows = size;
    mat->cols = size;

    return mat;
}

Matrix *random(int rows, int cols) {
    
    // input validation
    if (rows <= 0 || cols <= 0) { lml_error(__func__, "invalid dimensions"); return NULL; }

    // set generator seed
    srand((unsigned int)(time(NULL) * clock()));

    // allocate memory for matrix
    Matrix *mat = zeros(rows, cols);
    
    // initialize matrix elements to random values
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) { 
            mat->data[i][j] = (double)rand() / RAND_MAX; 
        }
    }

    return mat;
}

Matrix *matrix_from_array(int rows, int cols, double array[rows][cols]) {
    
    // input validation
    if (rows <= 0 || cols <= 0) { lml_error(__func__, "invalid dimensions"); return NULL; }
    if (array == NULL) { lml_error(__func__, "invalid array"); return NULL; }

    // allocate memory for matrix
    Matrix *mat = zeros(rows, cols);

    // copy the array data to the matrix
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) { 
            mat->data[i][j] = array[i][j]; 
        }
    }

    return mat;
}

// --------------- Retrieving Data ---------------

Matrix *copy(Matrix *mat) {
    
    // input validation
    if (mat == NULL) { lml_error(__func__, "invalid input matrix"); return NULL; }

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
    
    // input validation
    if (mat == NULL) { lml_error(__func__, "invalid input matrix"); return NULL; }
    if (row < 0 || row >= mat->rows) { lml_error(__func__, "invalid row index"); return NULL; }

    // create a new matrix to represent the row
    Matrix *row_matrix = zeros(1, mat->cols);

    // copy the row data into the new matrix
    for (int j = 0; j < mat->cols; j++) {
        row_matrix->data[0][j] = mat->data[row][j];
    }

    return row_matrix;
}

Matrix *get_col(Matrix *mat, int col) {
    
    // input validation
    if (mat == NULL) { lml_error(__func__, "invalid input matrix"); return NULL; }
    if (col < 0 || col >= mat->cols) { lml_error(__func__, "invalid column index"); return NULL; }

    // create a new matrix to represent the column
    Matrix *col_matrix = zeros(mat->rows, 1);

    // copy the column data into the new matrix
    for (int i = 0; i < mat->rows; i++) {
        col_matrix->data[i][0] = mat->data[i][col];
    }

    return col_matrix;
}

Matrix *get_lower(Matrix *mat) {
    
    // input validation
    if (mat == NULL) { lml_error(__func__, "invalid input matrix"); return NULL; }

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
    
    // input validation
    if (mat == NULL) { lml_error(__func__, "invalid input matrix"); return NULL; }

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
    
    // input validation
    if (mat == NULL) { lml_error(__func__, "invalid input matrix"); return NULL; }
    if (row < 0 || col < 0) { lml_error(__func__, "invalid indices"); return NULL; }
    if (rows <= 0 || cols <= 0) { lml_error(__func__, "invalid dimensions"); return NULL; }

    // create a new matrix for the submatrix
    Matrix *submatrix = zeros(rows, cols);

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
    
    // input validation
    if (mat == NULL) { lml_error(__func__, "invalid input matrix"); return 0.0; }
    if (mat->rows != mat->cols) { lml_error(__func__, "matrix is not NxN"); return 0.0; }

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
    
    // input validation
    if (mat == NULL) { lml_error(__func__, "invalid input matrix"); return NULL; }

    // create a new matrix with dimensions swapped
    Matrix *transposed = zeros(mat->cols, mat->rows);

    // copy elements from the original matrix to the transposed matrix
    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; j < mat->cols; j++) {
            transposed->data[j][i] = mat->data[i][j];
        }
    }

    return transposed;
}

Matrix *add(Matrix *mat1, Matrix *mat2) {
    
    // input validation
    if (mat1 == NULL || mat2 == NULL) { lml_error(__func__, "invalid input matrices"); return NULL; }
    if (mat1->rows != mat2->rows || mat1->cols != mat2->cols) { lml_error(__func__, "dimension mismatch"); return NULL; }

    // add matrices
    Matrix *result = zeros(mat1->rows, mat1->cols);
    for (int i = 0; i < result->rows; i++) {
        for (int j = 0; j < result->cols; j++) {
            result->data[i][j] = mat1->data[i][j] + mat2->data[i][j];
        }
    }

    return result;
}

Matrix *multiply(Matrix *mat1, Matrix *multiplier) {

    // input validation
    if (mat1 == NULL || multiplier == NULL) { lml_error(__func__, "invalid input matrices"); return NULL; }
    if (mat1->cols != multiplier->rows) { lml_error(__func__, "dimension mismatch"); return NULL; }

    // multiply matrices together
    Matrix *result = zeros(mat1->rows, multiplier->cols);
    for (int i = 0; i < mat1->rows; i++) {
        for (int j = 0; j < multiplier->cols; j++) {
            for (int k = 0; k < mat1->cols; k++) {
                result->data[i][j] += mat1->data[i][k] * multiplier->data[k][j];
            }
        }
    }

    return result;
}

Matrix *scalar_multiply(Matrix *mat, double scalar) {
    
    // input validation
    if (mat == NULL) { lml_error(__func__, "invalid input matrix"); return NULL; }

    // create a new matrix to store the result
    Matrix *result = copy(mat);

    // scale matrix
    for (int i = 0; i < result->rows; i++) {
        for (int j = 0; j < result->cols; j++) {
            result->data[i][j] *= scalar;
        }   
    }

    return result;
}

Matrix *scalar_add(Matrix *mat, double scalar) {
    
    // input validation
    if (mat == NULL) { lml_error(__func__, "invalid input matrix"); return NULL; }

    // create a new matrix to store the result
    Matrix *result = copy(mat);

    // scale matrix
    for (int i = 0; i < result->rows; i++) {
        for (int j = 0; j < result->cols; j++) {
            result->data[i][j] += scalar;
        }   
    }

    return result;
}

void LU_decompose(Matrix *mat, Matrix **L, Matrix **U) {
    
    // input validation
    if (mat == NULL) { lml_error(__func__, "invalid input matrix"); return; }    
    if (L == NULL || U == NULL) { lml_error(__func__, "invalid output matrices"); return; }    
    if (mat->rows != mat->cols) { lml_error(__func__, "matrix is not NxN"); return; }

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
    
    // input validation
    if (mat1 == NULL || mat2 == NULL) { lml_error(__func__, "invalid input matrices"); return NULL; }
    if (mat1->rows != mat1->cols) { lml_error(__func__, "matrix is not NxN"); return NULL; }
    if (mat1->rows != mat2->rows) { lml_error(__func__, "dimension mismatch"); return NULL; }

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
    
    // input validation
    if (mat == NULL) { lml_error(__func__, "invalid input matrix"); return NULL; }
    if (mat->rows != mat->cols) { lml_error(__func__, "matrix is not NxN"); return NULL; }

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

void map(Matrix *mat, double (*function)(double)) {
    
    // input validation
    if (mat == NULL) { lml_error(__func__, "invalid input matrix"); return; }
    if (function == NULL) { lml_error(__func__, "invalid input function"); return; }

    // apply the function element-wise to the matrix
    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; j < mat->cols; j++) {
            mat->data[i][j] = function(mat->data[i][j]);
        }
    }
}

void set_row(Matrix *mat, int row_index, Matrix *row_values) {
    
    // input validation
    if (mat == NULL || row_values == NULL) { lml_error(__func__, "invalid input matrix or row values"); return; }
    if (row_index < 0 || row_index >= mat->rows) { lml_error(__func__, "invalid row index"); return; }
    if (row_values->cols != mat->cols) { lml_error(__func__, "dimension mismatch"); return; }

    // copy the values from row_values to the matrix
    for (int j = 0; j < mat->cols; j++) {
        mat->data[row_index][j] = row_values->data[0][j];
    }
}

void set_col(Matrix *mat, int col_index, Matrix *col_values) {
    
    // input validation
    if (mat == NULL || col_values == NULL) { lml_error(__func__, "invalid input matrix or column values"); return; }
    if (col_index < 0 || col_index >= mat->cols) { lml_error(__func__, "invalid col index"); return; }
    if (col_values->rows != mat->rows) { lml_error(__func__, "dimension mismatch"); return; }

    // copy the values from col_values to the matrix
    for (int i = 0; i < mat->rows; i++) {
        mat->data[i][col_index] = col_values->data[i][0];
    }
}

void set_submatrix(Matrix *mat, int row, int col, Matrix *sub) {
    
    // input validation
    if (mat == NULL || sub == NULL) { lml_error(__func__, "invalid input matrix or submatrix"); return; }
    if (row < 0 || col < 0 || row + sub->rows > mat->rows || col + sub->cols > mat->cols) {
        lml_error(__func__, "invalid input submatrix dimensions"); return;
    }

    // copy values from the submatrix to the matrix
    for (int i = 0; i < sub->rows; i++) {
        for (int j = 0; j < sub->cols; j++) {
            mat->data[row + i][col + j] = sub->data[i][j];
        }
    }
}

void remove_row(Matrix *mat, int row) {
    
    // input validation
    if (mat == NULL) { lml_error(__func__, "invalid input matrix"); return; }
    if (row < 0 || row >= mat->rows) { lml_error(__func__, "invalid row index"); return; }

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
    
    // input validation
    if (mat == NULL) { lml_error(__func__, "invalid input matrix"); return; }
    if (col < 0 || col >= mat->cols) { lml_error(__func__, "invalid column index"); return; }

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
    
    // input validation
    if (mat == NULL || row_values == NULL) { lml_error(__func__, "invalid input matrix or row values"); return; }
    if (row < 0 || row > mat->rows) { lml_error(__func__, "invalid row index"); return; }
    if (row_values->cols != mat->cols) { lml_error(__func__, "dimension mismatch"); return; }

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
    
    // input validation
    if (mat == NULL || col_values == NULL) { lml_error(__func__, "invalid input matrix or column values"); return; }
    if (col < 0 || col > mat->cols) { lml_error(__func__, "invalid input matrix or column values"); return; }
    if (col_values->rows != mat->rows) { lml_error(__func__, "dimension mismatch"); return; }

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
    
    // input validation
    if (mat1 == NULL || mat2 == NULL) { lml_error(__func__, "invalid input matrices"); return; }
    if (mat1->cols != mat2->cols) { lml_error(__func__, "dimension mismatch"); return; }
    
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
    
    // input validation
    if (mat1 == NULL || mat2 == NULL) { lml_error(__func__, "invalid input matrices"); return; }
    if (mat1->rows != mat2->rows) { lml_error(__func__, "dimension mismatch"); return; }
    
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

    // input validation
    if (mat == NULL) { lml_error(__func__, "invalid input matrix"); return; }
    
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