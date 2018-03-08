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

int compressed::decomposeMatrix( comp_r_mat* DS , comp_r_mat* LUS , comp_r_mat* AS ){
	/* decomposes matrix AS into its diagonal elements (stored in DS)
	 * and its lower upper form (stored in LUS)
	 * create all three structures in main, and pass in as function handles
	 */

	// create a copy of AS in LUS
	copyMatrix( LUS , AS );

	// DS is a single vector of doubles containing the diagonal elements
	DS->noofRows = AS->noofRows;
	DS->noofCols = AS->noofCols;
	DS->noofVars = DS->noofRows;
	int row = 0; 
	DS->row_p.push_back(row);
	for ( int i = 0 ; i < DS->noofRows ; i++ ){
		// extract the diagonal elements to DS
		
		double value = retrieveElement( AS , i , i );

		if(value != 0.0){
			row++;
			DS->value.push_back(value);
			DS->col_id.push_back(i);
			//cout << "value: " << value << endl;
			//cout << "col: " << i << endl;
		}
		DS->row_p.push_back(row);
		
		//cout << "row: " << row << endl;

		// change the diagonal elements in LUS to 0
		changeElement( LUS , i , i , 0.0 );
	}

	// negate all the elements in LUS
	scalarMultiple( LUS , -1.0 );

	return 0;

}

int compressed::matrixProduct( comp_r_mat* result , comp_r_mat* A , comp_r_mat* vec ){
	int row = 0;
	result->row_p.push_back(row);
	for ( int i = 0 ; i < A->noofRows ; i++ ){
		double pdt = 0.0;
		for ( int j = 0 ; j < A->noofCols ; j++ ){
			pdt += retrieveElement( A , i , j )*retrieveElement( vec , j , j );
		}
		//cout << "pdt: " << pdt << endl;
		if ( pdt != 0.0 ){
			row++;
			result->value.push_back(pdt);
			result->col_id.push_back(0);
		}
		result->row_p.push_back(row);
	}
	return 0;
}

int compressed::jacobiSolver( comp_r_mat* X , comp_r_mat* DS , comp_r_mat* LUS , comp_r_mat* B ){
	int size = DS->noofRows;
	comp_r_mat matPdt;
	matPdt.noofRows = size;
	matPdt.noofCols = 1;
	matPdt.noofVars = size;
	//for ( int i = 0 ; i < size ; i ++ )	matPdt.value.push_back(0.0);

	matrixProduct(&matPdt , LUS , X );

	cout << endl << "X values: " << endl;
	for ( int i = 0 ; i < X->noofRows ; i++ ){
		cout << compressed::retrieveElement(X , i , 0 ) << "   ";
	}
	cout << endl;

	cout << endl << "product values: " << endl;

	for ( int i = 0 ; i < matPdt.noofRows ; i++ ){
		cout << compressed::retrieveElement(&matPdt , i , 0 ) << "   ";
	}
	cout << endl;

	comp_r_mat iter_result;
	iter_result.noofRows = X->noofRows;
	iter_result.noofCols = X->noofCols;
	iter_result.noofVars = X->noofVars;

	int row = 0;
	iter_result.row_p.push_back(row);
	cout << endl << "temp : ";
	for ( int i = 0 ; i < size ; i++ ){

		double dInv = 1.0/( retrieveElement(DS , i , i) );
		double temp = ( retrieveElement(B,i,0) + retrieveElement(&matPdt,i,0) ) * dInv;
		cout << temp << endl;
		if ( temp != 0.0 ){
			row++;
			iter_result.value.push_back(temp);
			iter_result.col_id.push_back(0);
		}
		iter_result.row_p.push_back(row);
	}
	/*cout << endl << "iter values: " << endl;
	for ( int i = 0 ; i < iter_result.noofRows ; i++ ){
		cout << compressed::retrieveElement(&iter_result , i , 0 ) << "   ";
	}
	cout << endl;
	cout << endl;*/
	(*X) = iter_result;

	return 0;

}

int compressed::calculateNorm( double& norm , vector<double>* v , vector<double>* Ax ){

}