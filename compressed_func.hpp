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
    	int noofRows;
    	int noofCols;
    	int noofVars;
	} comp_r_mat;

	typedef struct Compressed_Diagonal_Matrix{
		vector<double> value;
		int rank;		
	} comp_diag;

	void rowScale( comp_r_mat* A , int i , int j , double a );
	void rowPermute(comp_r_mat* A, int i, int j);
	double retrieveElement( comp_r_mat* input, int row_id, int col_id);
	comp_r_mat construct_compressed_matrix( vector<vector<double>>* input );
	comp_r_mat construct_compressed_matrix(vector<int>* i, vector<int>* j, vector<double>* val, int rowRank, int colRank);
	double productAx( comp_r_mat* A, vector<double>* x, vector<double>* b );
	void print_comp_r_mat( comp_r_mat* mat_a );
	bool check_sum( comp_r_mat* mat, vector<double>* vec );
	void reorderMat( comp_r_mat* input, comp_r_mat* reorder_A, comp_r_mat* reorder_B, int R, int C);
	void columnPermute(comp_r_mat* A, int col1, int col2);

	int changeElement( comp_r_mat* AS , int rowInd , int colInd , double newValue );
	int scalarMultiple( comp_r_mat* A , double scale );
	int copyMatrix( comp_r_mat* CS , comp_r_mat* AS );
	int decomposeMatrix( comp_diag* DS , comp_r_mat* LUS , comp_r_mat* AS );
	int matrixProduct ( comp_diag* result , comp_r_mat* A , comp_diag* vec);
	int jacobiSolver( comp_diag* X , comp_diag* DS , comp_r_mat* LUS , comp_diag* B );
	int calculateNorm( double& norm , vector<double>* v , vector<double>* Ax );

}
