// matrixmagic.h

#ifndef MATRIXMAGIC_H
#define MATRIXMAGIC_H

typedef struct {
    int rows;
    int cols;
    double **data;
} Matrix;

// matrix creation
Matrix *zeros(int rows, int cols);
Matrix *ones(int rows, int cols);
Matrix *identity(int size);
Matrix *random(int rows, int cols);
Matrix *matrix_from_array(int rows, int cols, double array[rows][cols]);

// data retrieval
double *get_row(Matrix *mat, int row);
double *get_col(Matrix *mat, int col);
double **as_array(Matrix *mat);
Matrix *copy(Matrix *mat);
Matrix *get_lower(Matrix *mat);
Matrix *get_upper(Matrix *mat);

// matrix operations
Matrix *transpose(Matrix *mat);
Matrix *inverse(Matrix *mat);
Matrix *multiply(Matrix *mat1, Matrix *mat2);
Matrix *solve_system(Matrix *mat1, Matrix *mat2);

// matrix editing
Matrix *map(Matrix *mat, double (*function)(double));
void set_row(Matrix *mat, int row_index, double *row_values);
void set_col(Matrix *mat, int col_index, double *col_values);
Matrix *remove_rows(Matrix *mat, int start_row, int num_rows);
Matrix *remove_cols(Matrix *mat, int start_col, int num_cols);
Matrix *insert_row(Matrix *mat, int row_index, double *row_values);
Matrix *insert_col(Matrix *mat, int col_index, double *col_values);

// displaying matrix
void display(Matrix *mat);

// deleting matrix
void destroy(Matrix *mat);

#endif /* MATRIXMAGIC_H */
