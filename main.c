/*
 * Example Usage of Lightweight Matrix Library (LML)
 * File: main.c
 * Author: James Bray
 *
 * This example demonstrates the usage of the Lightweight Matrix Library (LML)
 * by performing various operations on matrices such as generating matrices,
 * retrieving data, performing matrix operations, editing matrices, and displaying
 * results. The example showcases the functionality of the library and serves as
 * a reference for developers interested in using LML in their projects.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lml.h"

int main() {
    double data[3][3] = {
        {1, 2, 4},
        {3, 8, 14},
        {2, 6, 13} 
    };

    Matrix *mat1, *mat2, *mat3;
    
    printf("\n--------------- Generating Matrices ---------------\n");

    mat1 = zeros(2, 4); printf("\nGenerate zeros:\n"); display(mat1);
    mat1 = ones(3, 5); printf("\nGenerate ones:\n"); display(mat1);
    mat1 = identity(4); printf("\nGenerate identity:\n"); display(mat1);
    mat1 = random(6, 2); printf("\nGenerate random:\n"); display(mat1);
    mat1 = matrix_from_array(3, 3, data); printf("\nGenerate from array:\n"); display(mat1);

    printf("\n--------------- Retrieving Data ---------------\n");

    mat2 = get_row(mat1, 1); printf("\nGet row:\n"); display(mat2);
    mat2 = get_col(mat1, 1); printf("\nGet col:\n"); display(mat2);
    mat2 = copy(mat1); printf("\nCopying data:\n"); display(mat2);
    mat2 = get_lower(mat1); printf("\nGet lower:\n"); display(mat2);
    mat2 = get_upper(mat1); printf("\nGet upper:\n"); display(mat2);
    mat2 = get_submatrix(mat1, 1, 1, 2, 2); printf("\nGet submatrix\n"); display(mat2);

    printf("\n--------------- Matrix Operations ---------------\n");

    mat2 = transpose(mat1); printf("\nTranspose:\n"); display(mat2);

    Matrix *L; Matrix *U;
    LU_decompose(mat1, &L, &U);
    printf("\nLower:\n"); display(L);
    printf("\nUpper:\n"); display(U);
    release(L); release(U);

    printf("\nDeterminant: %f\n", det(mat1));

    printf("\n--------------- Matrix Editing ---------------\n");

    mat1 = identity(3); mat2 = ones(3, 3);
    scale(mat1, 2); printf("\nScale:\n"); display(mat1);
    shift(mat1, 0.5); printf("\nShift:\n"); display(mat1);

    mat3 = multiply(mat1, mat2); printf("\nMultiply:\n"); display(mat3);
    mat3 = add(mat1, mat2); printf("\nAdd:\n"); display(mat3);

    mat1 = zeros(4, 4); 
    mat2 = ones(1, 4); set_row(mat1, 1, mat2); 
    printf("\nSet row:\n"); display(mat1);
    mat2 = ones(4, 1); set_col(mat1, 1, mat2); 
    printf("\nSet col:\n"); display(mat1);

    remove_row(mat1, 1); printf("\nRemove row:\n"); display(mat1);
    remove_col(mat1, 1); printf("\nRemove column:\n"); display(mat1);

    // mat1 = random(5, 5); display(mat1);
    // mat2 = zeros(5, 1); insert_col(&mat1, 1, mat2);
    // printf("\nInsert col:\n"); display(mat1);

    // mat1 = random(5, 5); display(mat1);
    // mat2 = zeros(1, 5); insert_row(&mat1, 1, mat2);
    // printf("\nInsert row:\n"); display(mat1);

    release(mat1); release(mat2); release(mat3);

    return 0;
}