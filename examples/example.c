#include <stdio.h>
#include "matrixmagic.h"

int main() {
    // Example usage of MatrixMagiC library
    Matrix *mat1 = create_matrix(2, 2);
    // Populate mat1
    Matrix *mat2 = create_matrix(2, 2);
    // Populate mat2
    Matrix *result = add_matrices(mat1, mat2);
    // Display result or perform other operations
    free_matrix(mat1);
    free_matrix(mat2);
    free_matrix(result);
    return 0;
}