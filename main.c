// example.c

#include <stdio.h>
#include <stdlib.h>
#include "matrixmagic.h"

int main() {
    Matrix *mat = zeros(2, 4); display(mat);

    mat = ones(3, 5); display(mat);
    mat = identity(2); display(mat);
    mat = random(6, 2); display(mat);

    double data[3][3] = 
    {
        {1.23, 4.56, 2.34},
        {2.34, 3.45, 1.23},
        {4.56, 1.23, 5.03}
    };

    mat = matrix_from_array(3, 3, data); display(mat);

    destroy(mat);
    return 0;
}