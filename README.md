# ECE4960 - Programming Assignment 2 Full & Compressed Matrice and Jacobi Solver

***************************************************************************
## Part 1 - Full Matrix Implementation
***************************************************************************

**Overview:**
Implementation of structures and functions for full matrix operations.

**Documentation** 
Full matrices have to be created as a vector of vectors using the C++ stdlib.
The following functions are declared under the `full::` namespace.

- `retrieveElement` : overloaded function that can returns (i,j)th element in matrix or i-th element in vector

Examples:
```
vector< double > vec = { 1 , 2 , 3 , 4 , 5 };
vector< vector<double> > mat = { {1 , 2} , { 3 , 4 } };
double number;
number = full::retrieveElement( &vec , 2 );  //returns 3.0
number = full::retrieveElement( &mat , 0 , 1 ); //returns 2.0
```

- `print_full_mat( mat_ptr )` : print out whole matrix to screen
- `print_full_vec( vec_ptr )` : print out vector to screen

- `changeElement( mat_ptr , i , j , newValue )` : change (i,j)th element to newValue

- `copyMatrix( copy_ptr , original_ptr )` : create a copy of matrix in original_ptr in copy_ptr

- `scalarMultiple( mat_ptr , value )` : multiply elements in matrix by value

- `rowPermute( mat_ptr , i , j )` : swap rows i and j in matrix

- `rowScale( mat_ptr , i , j , a )` : set row_j = row_j + a * row_i in matrix

- `productAx( result_ptr , mat_ptr , vec_ptr )` : calculates vector product of mat and vec, and stores it in result; assumes mat is square and of the same rank as vec

- `calculateNorm( norm , vec_a , vec_b )` : calculates || vec_a - vec_b || and stores value in norm

- `decomposeMatrix( D_ptr , LU_ptr , mat_ptr )` : decomposes matrix into diagonals (stored in D) and lower-upper form (stored in LU)

- `jacobiSolver( X_ptr , D_ptr , LU_ptr , B_ptr )` : iteratively call this to find the solution to ( Ax = B ); solution stored in X_ptr

***************************************************************************
## Part 2 - Compressed Matrix Implementation
***************************************************************************

**Overview:**
Implementation of structures and functions for compressed matrix operations.

**Documentation** 
Compressed matrices have to be created using the compressed::comp_r_mat structure.
The following functions are declared under the `compressed::` namespace.

- `construct_compressed_matrix( )` : constructs a compressed matrix, given a vector of vectors

- `retrieveElement` : function that can returns (i,j)th element in matrix

Example:
```
vector< vector<double> > mat = { {1 , 0 , 0 } , { 0 , 1 , 4 } , { 0 , 0 , 3 }};

compressed::comp_r_mat AC = compressed::construct_compressed_matrix( &mat );
double number;
number = compressed::retrieveElement( &mat , 1 , 2 );  //returns 4.0
``` 

- `print_comp_r_mat( mat_ptr )` : print out compressed matrix representation to screen

- `changeElement( mat_ptr , i , j , newValue )` : change existing non-zero element at (i,j)th position to newValue

- `copyMatrix( copy_ptr , original_ptr )` : create a copy of matrix in original_ptr in copy_ptr

- `scalarMultiple( mat_ptr , value )` : multiply elements in matrix by value

- `rowPermute( mat_ptr , i , j )` : swap rows i and j in matrix

- `rowScale( mat_ptr , i , j , a )` : set row_j = row_j + a * row_i in matrix

- `productAx( result_ptr , mat_ptr , vec_ptr )` : calculates vector product of mat and vec, and stores it in result; assumes mat is square and of the same rank as vec.

- `calculateNorm( norm , vec_a , vec_b )` : calculates || vec_a - vec_b || and stores value in norm

- `decomposeMatrix( D_ptr , LU_ptr , mat_ptr )` : decomposes matrix into diagonals (stored in D) and lower-upper form (stored in LU)

- `jacobiSolver( X_ptr , D_ptr , LU_ptr , B_ptr )` : iteratively call this to find the solution to ( Ax = B ); solution stored in X_ptr.