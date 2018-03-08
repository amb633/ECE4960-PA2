#include <iostream> 
#include <cstdlib>
#include <cmath>
#include <fstream>

using namespace std;

namespace sparse{
	typedef struct sparse_matrix {
		double* value_array;
		int* rowPtr;
		int* colIdx;
		int noofRows;
		int noofCols;
		int noofVars;		
	} sparse_matrix;

	typedef struct sparse_diagonal_matrix {
		double* value_array;
		int noofVars;
	} sparse_diagonal_matrix;

	int retrieveElement( double& element , sparse_matrix* AS , int rowInd , int colIdx );
	int changeElement( sparse_matrix* AS , int rowInd , int colInd , double newValue );
	int copyMatrix( sparse_matrix* CS , sparse_matrix* AS );
	int rowPermute( sparse_matrix* AS , int row1 , int row2 );
	int rowScale( sparse_matrix* AS , int row1 , int row2 , double a );
	int matrixProduct( double* result , sparse_matrix* AS , double* vec );
	int decomposeMatrix( sparse_diagonal_matrix* DS , sparse_matrix* LUS , sparse_matrix* AS );
	int jacobiSolver( double* X , sparse_diagonal_matrix* DS , sparse_matrix* L , double* B );
	int calculateNorm( double& norm , double* v , double* Ax );
}