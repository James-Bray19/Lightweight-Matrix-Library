// example.c

#include <stdio.h>
#include <stdlib.h>
#include "matrixmagic.h"

int main() {
    
    // --------------- Generating Matrices ---------------
    printf("GENERATING MATRICES\n");

    printf("\nGenerate zeros:\n"); display(zeros(2, 4));
    printf("\nGenerate ones:\n"); display(ones(3, 5));
    printf("\nGenerate identity:\n"); display(identity(4));
    printf("\nGenerate random:\n"); display(random(6, 2));

    double data[4][4] = 
    {
        {1.23, 4.56, 2.34, 6.39},
        {2.34, 3.45, 1.23, 7.24},
        {5.58, 2.61, 4.47, 1.12},
        {4.56, 1.23, 5.03, 1.09}
    };

    Matrix *mat = matrix_from_array(4, 4, data);
    printf("\nGenerate from array:\n"); display(mat);

    // --------------- Retrieving Data ---------------
    printf("\n\nRETRIEVING DATA\n");

    printf("\nGet row:\n"); display(get_row(mat, 1));
    printf("\nGet col:\n"); display(get_col(mat, 1));

    printf("\nCopying data:\n");
    Matrix *mat2 = copy(mat); display(mat2);

    printf("\nGet lower:\n"); display(get_lower(mat));
    printf("\nGet upper:\n"); display(get_upper(mat));

    printf("\nGet submatrix\n"); display(get_submatrix(mat, 1, 1, 2, 2));

    destroy(mat);

    return 0;
}