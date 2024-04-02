NOTES:
    matrix_from_array needs to be optimised:
        function currently passes the array, needs to be optimised to pass pointers instead.
            ideal prototype:    
                Matrix *matrix_from_array(double **array, int rows, int cols);
            used prototype:
                Matrix *matrix_from_array(int rows, int cols, double array[rows][cols]);
        alternative solution: convert pointer type beforehand.
                double data[3][3] = 
                {
                    {1.23, 4.56, 2.34},
                    {2.34, 3.45, 1.23},
                    {4.56, 1.23, 5.03}
                };

                // Convert data to double **
                double *data_ptr[3];
                for (int i = 0; i < 3; ++i) {
                    data_ptr[i] = data[i];
                }               

                mat = matrix_from_array(data_ptr, 3, 3); display(mat);  

        

