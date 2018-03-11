//#ifndef full_mat_func_hpp
//#define full_mat_func_hpp

#include <iostream>
#include <cmath>
#include <vector>
#include <cstdlib>

using namespace std;

namespace full{

	/* UTILITY FUNCTIONS */

	// returns element at i,j in matrix AF / element i in vector VF 
	double retrieveElement( vector<vector<double>>* AF, int rowInd, int colInd);
	double retrieveElement( vector<double>* VF , int rowInd );

	// print the matrix AF / vector VF
	void print_full_mat( vector< vector<double>>* AF );
	void print_full_vec( vector<double>* VF );

	// set element at i,j in AF to newValue
	int changeElement( vector< vector<double>>* AF , int rowInd , int colInd , double newValue );

	// create a copy of matrix AF and store it in CF
	int copyMatrix( vector< vector<double>>* CF , vector< vector<double>>* AF );

	// multiply all (non-zero) elements in AF with scale
	int scalarMultiple( vector< vector<double>>* AF , double scale );

	// swaps row i,j in matrix AF
	int rowPermute( vector< vector <double>>* AF , int i , int j );

	// returns row_j = row_j + a * row_i
	int rowScale(vector< vector<double>>* AF, int i, int j, double a );

	// calculates the matrix product of AF and VF and stores it in result
	int productAx( vector<double>* result , vector< vector<double>>* AF , vector<double>* VF );

	// calculate norm of v - AX
	int calculateNorm( double& norm , vector<double>* v , vector<double>* Ax );	
	
	

	/* SOLVER FUNCTIONS */
	
	// decompose AF matrix into diagonal elements (stored in DF) and non-diagonal elements (stored in LUF)
	int decomposeMatrix( vector<double>* DF , vector< vector<double>>* LUF , vector< vector<double>>* AF );
	
	// jacobi solver function
	int jacobiSolver( vector<double>* X , vector<double>* DF , vector< vector<double>>* LUF , vector<double>* BF );

}