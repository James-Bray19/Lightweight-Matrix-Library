// example.c

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

    Matrix *mat1, *mat2;
    
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

    printf("\nDeterminant of Mat1: %f\n", det(mat1));

    printf("\n--------------- Matrix Editing ---------------\n");

    mat2 = ones(3, 4); printf("\nMultiply:\n"); display(multiply(mat1, mat2));
    printf("\nScale:\n"); display(scale(mat1, 3.0));

    release(mat1); release(mat2);

    return 0;
}