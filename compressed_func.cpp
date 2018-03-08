#include "compressed_func.hpp"

	void compressed::rowScale( comp_r_mat* A , int i , int j , double a ){

	}

	void compressed::rowPermute(comp_r_mat* A, int i, int j){

	}
	double compressed::retrieveElement( comp_r_mat* input, int row_id, int col_id){

	}
	comp_r_mat compressed::construct_compressed_matrix( vector<vector<double>> input ){

	}
	comp_r_mat compressed::construct_compressed_matrix(vector<int>* i, vector<int>* j, vector<double>* val, int rowRank, int colRank){

	}
	double compressed::productAx( comp_r_mat* A, vector<double>* x, vector<double>* b ){

	}
	void compressed::print_comp_r_mat( comp_r_mat* mat_a ){

	}
	bool compressed::check_sum( comp_r_mat* mat, vector<double>* vec ){

	}
	void compressed::reorderMat( comp_r_mat* input, comp_r_mat* reorder_A, comp_r_mat* reorder_B, int R, int C){

	}
	void compressed::columnPermute(comp_r_mat* A, int col1, int col2){

	}

	int compressed::changeElement( comp_r_mat* AS , int rowInd , int colInd , double newValue ){

	}
	int compressed::copyMatrix( comp_r_mat* CS , comp_r_mat* AS ){

	}
	int compressed::decomposeMatrix( comp_diag* DS , comp_r_mat* LUS , comp_r_mat* AS ){

	}
	int compressed::jacobiSolver( vector<double>* X , comp_diag* DS , comp_r_mat* LUS , vector<double>* B ){

	}
	int compressed::calculateNorm( double& norm , vector<double>* v , vector<double>* Ax ){

	}
