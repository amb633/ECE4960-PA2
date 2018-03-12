#ifndef compressed_func_hpp
#define compressed_func_hpp


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

	// returns element at i,j in matrix input
	double retrieveElement( comp_r_mat* input, int row_id, int col_id);

	// create a copy of matrix A and store it in C
	int copyMatrix( comp_r_mat* C , comp_r_mat* A );

    /* contructs compressed matrix from a double vector structure */
	comp_r_mat construct_compressed_matrix( vector<vector<double>>* input );
    
    /* contructs compressed matrix from a row and column index vector and value vector (memplux.mtx format) */
	comp_r_mat construct_compressed_matrix(vector<int>* i, vector<int>* j, vector<double>* val, int rowRank, int colRank);
	
    /* calculates the prodct of Matrix A with vector x and the results is in vector b */
    /* returns 0 if product is successful */
    /* returns 1 if Matrix A and vector x dimensions are not compatible */
    /* returns 2 if product results in any inf or ninf values */
    /* returns 3 if product results in any nan values */
    double productAx( vector<double>* b , comp_r_mat* A, vector<double>* x );
	
	/* prints the compressed matrix information */
	void print_comp_r_mat( comp_r_mat* mat_a );

	// multiplies all the values in a matrix by scale
	int scalarMultiple( comp_r_mat* A , double scale );

	// changes the value of an existing non-zero element
	int changeElement( comp_r_mat* A , int rowInd , int colInd , double newValue );

	// swaps rows i and j
	void rowPermute(comp_r_mat* A, int i, int j);

	// returns row_j = row_j + a*row_i
	void rowScale( comp_r_mat* A , int i , int j , double a );

	// calculate norm of v - Ax 
	int calculateNorm( double& norm , vector<double>* v , vector<double>* Ax );

	// decompose AS matrix into diagonal elements (stored in DS) and non-diagonal elements (stored in LUS)
	int decomposeMatrix( vector<double>* DS , comp_r_mat* LUS , comp_r_mat* AS );

	// jacobi solver function
	int jacobiSolver( vector<double>* X , vector<double>* DS , comp_r_mat* LUS , vector<double>* B );

}


#endif /* compressed_func_hpp */
