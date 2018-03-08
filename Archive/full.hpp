#include <iostream> 
#include <cstdlib>
#include <cmath>
#include <fstream>

using namespace std;

namespace full{
	typedef struct full_matrix {
		double** value_array;
		double noofRows;
		double noofCols;
	} full_matrix;

	typedef struct full_diagonal_matrix {
		double* value_array;
		int noofVars;
	} full_diagonal_matrix;

	int retrieveElement( double& element , full_matrix* AF , int rowInd , int colInd );
	int copyMatrix( full_matrix* CF , full_matrix* AF );
	int rowPermute( full_matrix* AF , int row1 , int row2 );
	int rowScale( full_matrix* AF , int row1 , int row2 , double a );
	int matrixProduct( double* result , full_matrix* AF , double* vec );
	int decomposeMatrix( full_diagonal_matrix* DF , full_matrix* LUF , full_matrix* AF );
}