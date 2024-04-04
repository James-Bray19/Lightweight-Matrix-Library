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
 * 
 * For more detailed documentation, see the README in the following repo:
 *      https://github.com/James-Bray19/Lightweight-Matrix-Library
 * 
 * Any contributions or issues should also be posted there.
 */

#include <stdio.h>
#include <stdlib.h>
#include "lml.h"

int main() {

    Matrix *mat1, *mat2, *mat3;


    printf("\n--------------- Generating Matrices ---------------\n");

    double data[3][3] = { {  1,  2,  4 },
                          {  3,  8, 14 },
                          {  2,  6, 13 } };

    printf("\nGenerate zeros:\n");
    mat1 = zeros(2, 4); display(mat1);

    printf("\nGenerate ones:\n");
    mat1 = ones(3, 5); display(mat1);

    printf("\nGenerate identity:\n");
    mat1 = identity(4); display(mat1);

    printf("\nGenerate random:\n");
    mat1 = random(6, 2); display(mat1);

    printf("\nGenerate from array:\n"); 
    mat1 = matrix_from_array(3, 3, data); display(mat1);



    printf("\n--------------- Retrieving Data ---------------\n");

    printf("\nThe following section will use this matrix:\n");
    mat1 = random(4, 4); display(mat1);

    printf("\nGet row:\n");
    mat2 = get_row(mat1, 1); display(mat2);

    printf("\nGet col:\n");
    mat2 = get_col(mat1, 1); display(mat2);

    printf("\nCopy of matrix:\n");
    mat2 = copy(mat1); display(mat2);

    printf("\nGet lower triangle:\n");
    mat2 = get_lower(mat1); display(mat2);

    printf("\nGet upper triangle:\n");
    mat2 = get_upper(mat1); display(mat2);

    printf("\nGet submatrix:\n");
    mat2 = get_submatrix(mat1, 1, 1, 2, 2); display(mat2);



    printf("\n--------------- Matrix Operations ---------------\n");

    Matrix *coeffs, *consts, *L, *U, *trans, *inv, *sol;

    printf("\nCoeffiecient Matrix:\n");
    coeffs = random(6, 6); display(coeffs);

    printf("\nConstants:\n");
    consts = random(6, 1); display(consts);

    printf("\nDeterminant of coefficient matrix:\n");
    printf("%8.2f\n", det(coeffs));

    printf("\nTranspose coefficients:\n");
    trans = transposed(coeffs); display(trans);

    printf("\nLU Decomposition:\n"); 
    LU_decompose(coeffs, &L, &U);
    printf("\nLower:\n"); display(L); 
    printf("\nUpper:\n"); display(U);

    printf("\nGaussian Solution:\n");
    sol = solve(coeffs, consts); display(sol);

    printf("\nInverse:\n");
    inv = inverse(coeffs); display(inv);

    printf("\nInverse Solution:\n");
    sol = multiply(inv, consts); display(sol);



    printf("\n--------------- Matrix Editing ---------------\n");

    mat1 = identity(3);

    printf("\nScale:\n");
    scale(mat1, 2); display(mat1);

    printf("\nShift:\n");
    shift(mat1, 0.5); display(mat1);
    
    printf("\nMultiply:\n");
    mat1 = identity(3); scale(mat1, 5); mat2 = ones(3, 3);
    mat3 = multiply(mat1, mat2); display(mat3);

    printf("\nAdd:\n");
    mat1 = identity(3); scale(mat1, 3); mat2 = zeros(3, 3);
    mat3 = add(mat1, mat2); display(mat3);

    mat1 = zeros(4, 4); 

    printf("\nSet row:\n");
    mat2 = ones(1, 4);
    set_row(mat1, 1, mat2); display(mat1);

    printf("\nSet col:\n"); 
    mat2 = ones(4, 1); 
    set_col(mat1, 1, mat2); display(mat1);

    printf("\nRemove row:\n");
    remove_row(mat1, 1); display(mat1);

    printf("\nRemove column:\n");
    remove_col(mat1, 1); display(mat1);

    printf("\nInsert row:\n"); 
    mat1 = random(4, 4); display(mat1); mat2 = zeros(1, 4);
    insert_row(mat1, 2, mat2); display(mat1);

    printf("\nInsert col:\n"); 
    mat1 = random(4, 4); display(mat1); mat2 = zeros(4, 1);
    insert_col(mat1, 2, mat2); display(mat1);

    printf("\nAppend rows:\n"); 
    mat1 = random(3, 3); display(mat1); mat2 = zeros(2, 3);
    append_rows(mat1, mat2); display(mat1);

    printf("\nAppend cols:\n"); 
    mat1 = random(3, 3); display(mat1); mat2 = zeros(3, 2);
    append_cols(mat1, mat2); display(mat1);

    printf("\n--------------- Program Finished ---------------\n");

    // release memory after processing
    release(mat1); release(mat2); release(mat3);
    release(coeffs); release(consts); release(L); release(U);
    release(trans); release(inv); release(sol);

    return 0;
}