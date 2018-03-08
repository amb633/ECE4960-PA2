#include "compressed_func.hpp"

void compressed::rowScale( comp_r_mat* A , int i , int j , double a ){

}

void compressed::rowPermute(comp_r_mat* A, int i, int j){

}

double compressed::retrieveElement( comp_r_mat* input, int row_id, int col_id){
	double element = 0;
    int row_non_zero_start = input->row_p[row_id];
    int row_non_zero_end = input->row_p[row_id + 1];
    for(int p = row_non_zero_start; p< row_non_zero_end; p++){
        if( input->col_id[p] == col_id ){
            element = input->value[p];
            break;
        }
    }
    return element;
}

compressed::comp_r_mat compressed::construct_compressed_matrix( vector<vector<double>>* input ){
	comp_r_mat mat_A;
	int rows = (*input).size();
	int columns = (*input)[0].size();
	int non_zeros = 0;
	mat_A.row_p.push_back(non_zeros);

	for ( int i = 0 ; i < rows ; i++ ){
		for ( int j = 0 ; j < columns ; j++ ){
			if ((*input)[i][j] !=0){
				mat_A.value.push_back((*input)[i][j]);
				mat_A.col_id.push_back(j);
				non_zeros++;
			}
		}
		mat_A.row_p.push_back(non_zeros);
	}
	return mat_A;
}

compressed::comp_r_mat compressed::construct_compressed_matrix(vector<int>* i, vector<int>* j, vector<double>* val, int rowRank, int colRank){

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

int compressed::changeElement( comp_r_mat* A , int rowInd , int colInd , double newValue ){
	/* changes the value of an existing non-zero element
	 * returns 0 if the element was successfully changed 
	 * returns -1 if the element is currently a zero, and couldn't be changed 
	 */
	int start = A->row_p[rowInd];
	int end = A->row_p[rowInd+1];
	for ( int i = start ; i < end ; i++ ){
		if ( A->col_id[i] == colInd ){
			A->value[i] = newValue;
			return 0;
		}
	}
	return -1;
}

int compressed::scalarMultiple( comp_r_mat* A , double scale ){
	/*  scalar operations (e.g. want to multiply all elements in matrix by 2);
	 *  so the row and column order are not important, just muliply everything
	 */
	for ( int i = 0 ; i < A->noofVars ; i++ ){
		if( A->value[i] != 0 )	A->value[i] = scale*A->value[i];
	}
	return 0;
}

int compressed::copyMatrix( comp_r_mat* C , comp_r_mat* A ){
	/* creates a copy of matrix A in matrix C
	 * create a pointer to both in the main, and pass the pointers as function handles
	 * need to do a deep copy if using dynamic memory allocation
	 */
	C->noofRows = A->noofRows;
	C->noofCols = A->noofCols;
	C->noofVars = A->noofVars;
	C->value = A->value;
	C->row_p = A->row_p;
	C->col_id = A->col_id;
	return 0;

}

int compressed::decomposeMatrix( comp_diag* DS , comp_r_mat* LUS , comp_r_mat* AS ){
	/* decomposes matrix AS into its diagonal elements (stored in DS)
	 * and its lower upper form (stored in LUS)
	 * create all three structures in main, and pass in as function handles
	 */

	// create a copy of AS in LUS
	copyMatrix( LUS , AS );

	// DS is a single vector of doubles containing the diagonal elements
	DS->rank = AS->noofRows;

	for ( int i = 0 ; i < DS->rank ; i++ ){
		// extract the diagonal elements to DS
		DS->value.push_back(retrieveElement( AS , i , i ));
		// change the diagonal elements in LUS to 0
		changeElement( LUS , i , i , 0.0 );
	}

	// negate all the elements in LUS
	scalarMultiple( LUS , -1.0 );

	return 0;

}

int compressed::matrixProduct( comp_diag* result , comp_r_mat* A , comp_diag* vec ){
	for ( int i = 0 ; i < A->noofRows ; i++ ){
		for ( int j = 0 ; j < A->noofCols ; j++ ){
			double temp = retrieveElement( A , i , j );
			result->value[i] += temp*(vec->value[j]);
			//cout << result->value[i] << endl;
		}
	}
	return 0;
}

int compressed::jacobiSolver( comp_diag* X , comp_diag* DS , comp_r_mat* LUS , comp_diag* B ){
	int size = DS->rank;
	comp_diag matPdt;
	matPdt.rank = size;
	for ( int i = 0 ; i < size ; i ++ )	matPdt.value.push_back(0.0);

	matrixProduct(&matPdt , LUS , X );
	
	for ( int i = 0 ; i < size ; i++ ){
		double dInv = 1.0/DS->value[i];
		X->value[i] = (B->value[i] + matPdt.value[i])*dInv;
	}

	return 0;

}

int compressed::calculateNorm( double& norm , vector<double>* v , vector<double>* Ax ){

}