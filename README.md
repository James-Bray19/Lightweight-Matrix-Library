# Lightweight Matrix Library (LML)

## Overview

**Lightweight Matrix Library (LML)** is a collection of functions designed for efficient manipulation and operation on matrices. This library is specifically tailored for embedded systems, offering functionality for generating matrices, retrieving data, performing matrix operations, editing matrices, and various miscellaneous functions.

Please see `docs/documentation.md` for a more detailed explanation of this library's functionality.

## Usage

To use the Lightweight Matrix Library in your project, follow these steps:

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

## Contributions

Any contributions or issues related to the Lightweight Matrix Library can be posted in the issues section. Your feedback is greatly appreciated!


