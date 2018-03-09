//#ifndef full_mat_func_hpp
//#define full_mat_func_hpp

#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <stdio.h>

using namespace std;

namespace full{

	/*typedef struct Full_Matrix{
		vector<double> value;
    	vector<size_t> row_p;
    	vector<int> col_id;
	} full_mat;

	typedef struct Full_Diagonal_Matrix{
		vector<double> value;
		vector<size_t> rank;		
	} full_diag;*/

	double retrieveElement( vector<vector<double>>* input, int row_id, int col_id);
	double productAx( vector< vector<double>>* input, vector<double>* x, vector<double>* b);
	void rowScale(vector< vector<double>>* input, int i, int j, double a );
	void rowPermute(vector< vector<double>>* input, int i, int j);

	void print_full_mat( vector< vector<double>>* input );
	//bool check_sum( vector< vector<double>>* mat, vector<double>* vec );
	void reorderMat( vector< vector<double>>* input, vector< vector<double>>* reorder_A, vector< vector<double>>* reorder_B, int R, int C);
	void columnPermute(vector< vector<double>>* A, int col1, int col2);

	int changeElement( vector< vector<double>>* A , int rowInd , int colInd , double newValue );
	int scalarMultiple( vector< vector<double>>* A , double scale );
	int copyMatrix( vector< vector<double>>* C , vector< vector<double>>* A );
	int decomposeMatrix( vector<double>* D , vector< vector<double>>* LU , vector< vector<double>>* A );
	int jacobiSolver( vector<double>* X , vector<double>* D , vector< vector<double>>* LU , vector<double>* B );
	int calculateNorm( double& norm , vector<double>* v , vector<double>* Ax );

}