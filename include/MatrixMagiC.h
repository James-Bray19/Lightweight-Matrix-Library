// matrixmagic.h

#ifndef MATRIXMAGIC_H
#define MATRIXMAGIC_H

typedef struct {
    int rows;
    int cols;
    double **data;
} Matrix;

Matrix *zeros(int rows, int cols);
Matrix *ones(int rows, int cols);
Matrix *identity(int size);
Matrix *random(int rows, int cols);
Matrix *matrix_from_array(int rows, int cols, double array[rows][cols]);

Matrix *get_row(Matrix *mat, int row);
Matrix *get_col(Matrix *mat, int col);
Matrix *copy(Matrix *mat);
Matrix *get_lower(Matrix *mat);
Matrix *get_upper(Matrix *mat);
Matrix *get_submatrix(Matrix *mat, int row, int col, int rows, int cols);

Matrix *transpose(Matrix *mat);
Matrix *inverse(Matrix *mat);
Matrix *multiply(Matrix *mat1, Matrix *mat2);
Matrix *solve_system(Matrix *mat1, Matrix *mat2);

Matrix *map(Matrix *mat, double (*function)(double));
void set_row(Matrix *mat, int row_index, double *row_values);
void set_col(Matrix *mat, int col_index, double *col_values);
Matrix *remove_rows(Matrix *mat, int start_row, int num_rows);
Matrix *remove_cols(Matrix *mat, int start_col, int num_cols);
Matrix *insert_row(Matrix *mat, int row_index, double *row_values);
Matrix *insert_col(Matrix *mat, int col_index, double *col_values);

void display(Matrix *mat);
void destroy(Matrix *mat);

#endif /* MATRIXMAGIC_H */