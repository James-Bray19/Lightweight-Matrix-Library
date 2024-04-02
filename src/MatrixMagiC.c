#include "MatrixMagiC.h"
#include <stdio.h>
#include <stdlib.h>

Matrix *createMatrix(int rows, int cols) {
    Matrix *matrix = (Matrix *)malloc(sizeof(Matrix));
    if (matrix == NULL) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    
    matrix->rows = rows;
    matrix->cols = cols;
    matrix->data = (double **)malloc(rows * sizeof(double *));
    if (matrix->data == NULL) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < rows; i++) {
        matrix->data[i] = (double *)malloc(cols * sizeof(double));
        if (matrix->data[i] == NULL) {
            fprintf(stderr, "Error: Memory allocation failed\n");
            exit(EXIT_FAILURE);
        }
    }
    
    // Initialize matrix elements to zero
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            matrix->data[i][j] = 0.0;
        }
    }
    
    return matrix;
}

void destroyMatrix(Matrix *matrix) {
    if (matrix != NULL) {
        for (int i = 0; i < matrix->rows; i++) {
            free(matrix->data[i]);
        }
        free(matrix->data);
        free(matrix);
    }
}

void printMatrix(Matrix *matrix) {
    if (matrix == NULL) {
        printf("Matrix is NULL\n");
        return;
    }
    printf("Matrix (%d x %d):\n", matrix->rows, matrix->cols);
    for (int i = 0; i < matrix->rows; i++) {
        for (int j = 0; j < matrix->cols; j++) {
            printf("%.2f ", matrix->data[i][j]);
        }
        printf("\n");
    }
}