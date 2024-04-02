// matrixmagic.h

#ifndef MATRIXMAGIC_H
#define MATRIXMAGIC_H

typedef struct {
    int rows;
    int cols;
    double **data;
} Matrix;


/*
// function to create a matrix initialized with zeros
// parameters:
// - rows: number of rows in the matrix
// - cols: number of columns in the matrix
// returns:
// - pointer to the created matrix
*/
Matrix *zeros(int rows, int cols);

/*
// function to create a matrix initialized with ones
// parameters:
// - rows: number of rows in the matrix
// - cols: number of columns in the matrix
// returns:
// - pointer to the created matrix
*/
Matrix *ones(int rows, int cols);

/*
// function to create an identity matrix
// parameters:
// - size: size of the square identity matrix (number of rows/columns)
// returns:
// - pointer to the created identity matrix
*/
Matrix *identity(int size);

/*
// function to create a matrix initialized with random values
// parameters:
// - rows: number of rows in the matrix
// - cols: number of columns in the matrix
// returns:
// - pointer to the created matrix
*/
Matrix *random(int rows, int cols);

/*
// function to create a matrix from a 2D array of doubles
// parameters:
// - rows: number of rows in the matrix
// - cols: number of columns in the matrix
// - array: 2D array of doubles to populate the matrix
// returns:
// - pointer to the created matrix
*/
Matrix *matrix_from_array(int rows, int cols, double array[rows][cols]);


/*
// function to get a pointer to a specific row of a matrix
// parameters:
// - mat: pointer to the matrix
// - row: index of the row to retrieve (0-based indexing)
// returns:
// - pointer to the matrix representing the requested row
*/
Matrix *get_row(Matrix *mat, int row);

/*
// function to get a pointer to a specific column of a matrix
// parameters:
// - mat: pointer to the matrix
// - col: index of the column to retrieve (0-based indexing)
// returns:
// - pointer to the matrix representing the requested column
*/
Matrix *get_col(Matrix *mat, int col);

/*
// function to create a copy of a matrix
// parameters:
// - mat: pointer to the matrix to be copied
// returns:
// - pointer to the copied matrix
*/
Matrix *copy(Matrix *mat);

/*
// function to extract the lower triangular part of a matrix
// parameters:
// - mat: pointer to the matrix
// returns:
// - pointer to the lower triangular part of the matrix
*/
Matrix *get_lower(Matrix *mat);

/*
// function to extract the upper triangular part of a matrix
// parameters:
// - mat: pointer to the matrix
// returns:
// - pointer to the upper triangular part of the matrix
*/
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