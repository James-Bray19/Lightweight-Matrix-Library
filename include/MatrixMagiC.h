#ifndef MATRIX_MAGIC_H
#define MATRIX_MAGIC_H

// Define a structure for matrices
typedef struct {
    int rows;
    int cols;
    double **data; // 2D array to store matrix elements
} Matrix;

// Function declarations for matrix operations
Matrix *createMatrix(int rows, int cols);
void destroyMatrix(Matrix *matrix);
void printMatrix(Matrix *matrix);
Matrix *add(Matrix *matrix1, Matrix *matrix2);
Matrix *appendAsColumn(Matrix *matrix1, Matrix *matrix2);
Matrix *appendAsRow(Matrix *matrix1, Matrix *matrix2);
// Add more function declarations for other matrix operations...

#endif /* MATRIX_MAGIC_H */
