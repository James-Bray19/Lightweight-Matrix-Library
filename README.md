# Lightweight Matrix Library (LML)

## Overview

**Lightweight Matrix Library (LML)** is a collection of functions designed for efficient manipulation and operation on matrices. This library is specifically tailored for embedded systems, offering functionality for generating matrices, retrieving data, performing matrix operations, editing matrices, and various miscellaneous functions.

## Usage

### Using the library

1. **Include Header File**: 
   - In your C program, include the `lml.h` header file to gain access to the library's functions and data structures.
   ```c
   #include "lml.h"
   ```

2. **Compile Your Program**: 
   - Use the provided Makefile to compile your program. If you want to make changes to the library code (`lml.c`), you can do so and then recompile your program using the Makefile.
   ```bash
   make
   ```

3. **Link Against the Library**: 
   - When compiling your program, ensure that you link against the library files (`lib/liblml.so`).
   ```bash
   gcc -o your_program your_program.c -Llib -llml
   ```

4. **Run Your Program**: 
   - After successfully compiling your program, you can run the executable as usual.
   ```bash
   ./my_program
   ```

### Modifying the library

1. **Modify the Source and Header Files**:
   - Make the necessary changes to the source files (`lml.c`) and header file (`lml.h`) according to your requirements.

2. **Compile the Library**:
   - Use the provided Makefile included in the repository to compile the library. Please edit the Makefile should you require a different library format.

3. **Test Your Changes**:
   - Before submitting your changes, it's essential to test them thoroughly to ensure that they work as expected.

4. **Submit a Pull Request**:
   - Once you are satisfied with your changes and have tested them thoroughly, submit a pull request to the repository.

## Matrix Structure
```c
typedef struct {
    int rows;
    int cols;
    double **data;
} Matrix;
```

The `Matrix` structure represents a two-dimensional matrix and consists of the following members:

- `rows`: An integer representing the number of rows in the matrix.
- `cols`: An integer representing the number of columns in the matrix.
- `data`: A pointer to a two-dimensional array of double values, which stores the actual data of the matrix.

The `Matrix` structure is used as the fundamental data type in the LML library to perform various matrix operations and manipulations. It provides a flexible and efficient way to work with matrices in C programming.

## Functions

#### `Matrix *zeros(int rows, int cols);`
> This function returns a matrix filled with zeros of the specified size.
> - `int rows`: Number of rows in the matrix.
> - `int cols`: Number of columns in the matrix.

#### `Matrix *ones(int rows, int cols);`
> This function returns a matrix filled with ones of the specified size.
> - `int rows`: Number of rows in the matrix.
> - `int cols`: Number of columns in the matrix.

#### `Matrix *identity(int size);`
> This function returns an identity matrix of the specified size.
> - `int size`: Size of the identity matrix.

#### `Matrix *random(int rows, int cols);`
> This function returns a matrix filled with random values of the specified size.
> - `int rows`: Number of rows in the matrix.
> - `int cols`: Number of columns in the matrix.

#### `Matrix *matrix_from_array(int rows, int cols, double array[rows][cols]);`
> This function creates a matrix from a 2D array.
> - `int rows`: Number of rows in the matrix.
> - `int cols`: Number of columns in the matrix.
> - `double array[rows][cols]`: Input array to create the matrix from.

#### `Matrix *copy(Matrix *mat);`
> This function creates a copy of the input matrix.
> - `Matrix *mat`: Input matrix to be copied.

#### `Matrix *get_row(Matrix *mat, int row);`
> This function returns the specified row of the matrix.
> - `Matrix *mat`: Input matrix.
> - `int row`: Row index.

#### `Matrix *get_col(Matrix *mat, int col);`
> This function returns the specified column of the matrix.
> - `Matrix *mat`: Input matrix.
> - `int col`: Column index.

#### `Matrix *get_lower(Matrix *mat);`
> This function returns the lower triangular matrix of the input matrix.
> - `Matrix *mat`: Input matrix.

#### `Matrix *get_upper(Matrix *mat);`
> This function returns the upper triangular matrix of the input matrix.
> - `Matrix *mat`: Input matrix.

#### `Matrix *get_submatrix(Matrix *mat, int row, int col, int rows, int cols);`
> This function returns a submatrix of the input matrix.
> - `Matrix *mat`: Input matrix.
> - `int row`: Starting row index.
> - `int col`: Starting column index.
> - `int rows`: Number of rows in the submatrix.
> - `int cols`: Number of columns in the submatrix.

#### `double det(Matrix *mat);`
> This function returns the determinant of the input matrix.
> - `Matrix *mat`: Input matrix.

#### `Matrix *transpose(Matrix *mat);`
> This function returns the transposed matrix.
> - `Matrix *mat`: Input matrix.

#### `Matrix *add(Matrix *mat1, Matrix *mat2);`
> This function returns the result of element-wise addition of two matrices.
> - `Matrix *mat1`: First matrix.
> - `Matrix *mat2`: Second matrix.

#### `Matrix *multiply(Matrix *mat1, Matrix *mat2);`
> This function returns the result of matrix multiplication of two matrices.
> - `Matrix *mat1`: First matrix.
> - `Matrix *mat2`: Second matrix.

#### `void *LU_decompose(Matrix *mat, Matrix **L, Matrix **U);`
> This function decomposes the input matrix into Lower and Upper triangular matrices.
> - `Matrix *mat`: Input matrix.
> - `Matrix **L`: Output Lower triangular matrix.
> - `Matrix **U`: Output Upper triangular matrix.

#### `Matrix *solve(Matrix *mat1, Matrix *mat2);`
> This function solves a system of linear equations represented by matrices.
> - `Matrix *mat1`: Coefficient matrix.
> - `Matrix *mat2`: Constant matrix.

#### `Matrix *inverse(Matrix *mat);`
> This function returns the inverse of the input matrix.
> - `Matrix *mat`: Input matrix.

#### `void scale(Matrix *mat, double scalar);`
> This function scales the matrix by a scalar value.
> - `Matrix *mat`: Input matrix.
> - `double scalar`: Scalar value.

#### `void shift(Matrix *mat, double scalar);`
> This function shifts the matrix by adding a scalar value to each element.
> - `Matrix *mat`: Input matrix.
> - `double scalar`: Scalar value.

#### `void map(Matrix *mat, double (*function)(double));`
> This function applies a function element-wise to the matrix.
> - `Matrix *mat`: Input matrix.
> - `double (*function)(double)`: Function pointer to apply.

#### `void set_row(Matrix *mat, int row_index, Matrix *row_values);`
> This function sets the values of a specific row in the matrix.
> - `Matrix *mat`: Input matrix.
> - `int row_index`: Row index.
> - `Matrix *row_values`: Matrix of row values.

#### `void set_col(Matrix *mat, int col_index, Matrix *col_values);`
> This function sets the values of a specific column in the matrix.
> - `Matrix *mat`: Input matrix.
> - `int col_index`: Column index.
> - `Matrix *col_values`: Matrix of column values.

#### `void remove_row(Matrix *mat, int row);`
> This function removes a row from the matrix.
> - `Matrix *mat`: Input matrix.
> - `int row`: Index of row to remove.

#### `void remove_col(Matrix *mat, int col);`
> This function removes a column from the matrix.
> - `Matrix *mat`: Input matrix.
> - `int col`: Index of column to remove.

#### `void insert_row(Matrix *mat, int row, Matrix *row_values);`
> This function inserts a new row with the provided values at the specified row index.
> - `Matrix *mat`: Destination matrix.
> - `int row`: Row index.
> - `Matrix *row_values`: Array of row values.

#### `void insert_col(Matrix *mat, int col, Matrix *col_values);`
> This function inserts a new column with the provided values at the specified column index.
> - `Matrix *mat`: Destination matrix.
> - `int col`: Column index.
> - `Matrix *col_values`: Array of column values.

#### `void append_rows(Matrix *mat1, Matrix *mat2);`
> This function appends rows from mat2 to mat1.
> - `Matrix *mat1`: Destination matrix.
> - `Matrix *mat2`: Source matrix.

#### `void append_cols(Matrix *mat1, Matrix *mat2);`
> This function appends columns from mat2 to mat1.
> - `Matrix *mat1`: Destination matrix.
> - `Matrix *mat2`: Source matrix.

#### `void display(Matrix *mat);`
> This function displays the content of the matrix.
> - `Matrix *mat`: Input matrix.

#### `void release(Matrix *mat);`
> This function releases memory allocated for the matrix.
> - `Matrix *mat`: Input matrix.

## Contributions

Any contributions or issues related to the Lightweight Matrix Library are welcome. Your feedback is greatly appreciated!


