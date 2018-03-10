//#ifndef compressed_mat_func_hpp
//#define compressed_mat_func_hpp


#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <fstream>

using namespace std;

namespace compressed{

	typedef struct Compressed_Row_Matrix{
		vector<double> value;
    	vector<int> row_p;
    	vector<int> col_id;
	} comp_r_mat;

	/* UTILITY FUNCTIONS */

	// creating the compressed matrices
	comp_r_mat construct_compressed_matrix( vector<vector<double>>* input );
	comp_r_mat construct_compressed_matrix(vector<int>* i, vector<int>* j, vector<double>* val, int rowRank, int colRank);

	// returns element at i,j in matrix A
	double retrieveElement( comp_r_mat* input, int row_id, int col_id);

	// create a copy of matrix A and store it in C
	int copyMatrix( comp_r_mat* C , comp_r_mat* A );

	// prints the compressed matrix information
	void print_comp_r_mat( comp_r_mat* mat_a );

	// multiplies all the values in a matrix by scale
	int scalarMultiple( comp_r_mat* A , double scale );

	// changes the value of an existing non-zero element */
	int changeElement( comp_r_mat* A , int rowInd , int colInd , double newValue );

	// swaps rows i and j
	void rowPermute(comp_r_mat* A, int i, int j);

	// returns row_j = row_j + a*row_i
	void rowScale( comp_r_mat* A , int i , int j , double a );
	
	// calculates dot product ( do we need this function? )
	double productAx( comp_r_mat* A, vector<double>* x, vector<double>* b );

	// calculates the matrix product of A and vec and stores it in result
	int matrixProduct ( vector<double>* result , comp_r_mat* A , vector<double>* vec);

	// calculate norm of v - Ax 
	int calculateNorm( double& norm , vector<double>* v , vector<double>* Ax );

	bool check_sum( comp_r_mat* mat, vector<double>* vec );
	void reorderMat( comp_r_mat* input, comp_r_mat* reorder_A, comp_r_mat* reorder_B, int R, int C);
	void columnPermute(comp_r_mat* A, int col1, int col2);	

	// decompose AS matrix into diagonal elements (stored in DS) and non-diagonal elements (stored in LUS)
	int decomposeMatrix( vector<double>* DS , comp_r_mat* LUS , comp_r_mat* AS );

	// jacobi solver function
	int jacobiSolver( vector<double>* X , vector<double>* DS , comp_r_mat* LUS , vector<double>* B );

}