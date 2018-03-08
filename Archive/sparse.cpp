#include "sparse.hpp"

int sparse::retrieveElement( double& element , sparse_matrix* matrix , int rowInd , int colInd ){
	/* find out how elements until previous row */
	int nnzCount = matrix->rowPtr[rowInd];
	/* find out how many (cumulative) elements in current row */
	int nnz = matrix->rowPtr[rowInd+1];

	for ( int i = nnzCount ; i < nnz ; i++ ){
		if ( matrix->colIdx[i] == colInd){
			element = matrix->value_array[i];
			return 0;
		}
	}
	element = 0.0;
	return 0;
}

int sparse::changeElement( sparse_matrix* matrix , int rowInd , int colInd , double newValue ){

	int nnzCount = matrix->rowPtr[rowInd];
	int nnz = matrix->rowPtr[rowInd+1];

	for ( int i = nnzCount ; i < nnz ; i++ ){
		if ( matrix->colIdx[i] == colInd ){
			matrix->value_array[i] = newValue;
			return 0;
		}
	}
	return 0;
}

int sparse::rowPermute( sparse_matrix* A , int row1 , int row2 ){
	
	/* ensure that row1 is the smaller one */
	if ( row2 < row1 ){
		int temp = row2;
		row2 = row1;
		row1 = temp;
	}

	/* Define some counters and limits */
	int nnCount_1 = A->rowPtr[row1];
	int nnLimit_1 = A->rowPtr[row1+1];
	int nnCount_2 = A->rowPtr[row2];
	int nnLimit_2 = A->rowPtr[row2+1];

	/* Create a new sparse matrix */
	sparse::sparse_matrix* temp_matrix = new sparse::sparse_matrix;
	temp_matrix->value_array = A->value_array;
	temp_matrix->colIdx = A->colIdx;
	temp_matrix->rowPtr = A->rowPtr;
	temp_matrix->noofRows = A->noofRows;
	temp_matrix->noofCols = A->noofCols;
	temp_matrix->noofVars = A->noofVars;

	int count = 0;
	int rowCount = 0;

	A->value_array = new double[temp_matrix->noofVars];
	A->colIdx = new int[temp_matrix->noofVars];
	A->rowPtr = new int[temp_matrix->noofRows+1];

	/* copy elements upto row 1 */
	for ( int i = temp_matrix->rowPtr[0] ; i < nnCount_1 ; i++ ) {
		A->value_array[count] = temp_matrix->value_array[i];
		A->colIdx[count] = temp_matrix->colIdx[i];
		count++;
	}

	for ( int i = 0 ; i <= row1 ; i++ ) {
		A->rowPtr[rowCount] = temp_matrix->rowPtr[i];
		rowCount++;
	}

	/* copy row 2 */
	for ( int i = nnCount_2 ; i < nnLimit_2 ; i++ ){
		A->value_array[count] = temp_matrix->value_array[i];
		A->colIdx[count] = temp_matrix->colIdx[i];
		count++;
	}

	A->rowPtr[rowCount] = A->rowPtr[rowCount-1] + (nnLimit_2 - nnCount_2);
	rowCount++;

	/* copy elements between row 1 and row 2 */
	for ( int i = nnLimit_1 ; i < nnCount_2 ; i++ ){
		A->value_array[count] = temp_matrix->value_array[i];
		A->colIdx[count] = temp_matrix->colIdx[i];
		count++;
	}

	for ( int i = row1+1 ; i < row2 ; i++ ) {
		A->rowPtr[rowCount] = A->rowPtr[rowCount-1]+(temp_matrix->rowPtr[i+1] - temp_matrix->rowPtr[i]);
		rowCount++;
	}

	/* copy row 1 */
	for ( int i = nnCount_1 ; i < nnLimit_1 ; i++ ) {
		A->value_array[count] = temp_matrix->value_array[i];
		A->colIdx[count] = temp_matrix->colIdx[i];
		count++;
	}

	A->rowPtr[rowCount] = A->rowPtr[rowCount-1] + (nnLimit_1 - nnCount_1);
	rowCount++;

	/* copy elements from row 2 to end */
	for ( int i = nnLimit_2 ; i < temp_matrix->rowPtr[A->noofRows] ; i++ ){
		A->value_array[count] = temp_matrix->value_array[i];
		A->colIdx[count] = temp_matrix->colIdx[i];
		count++;
	}

	for ( int i = row2+1 ; i < temp_matrix->noofRows ; i++ ) {
		A->rowPtr[rowCount] = A->rowPtr[rowCount-1]+(temp_matrix->rowPtr[i+1] - temp_matrix->rowPtr[i]);
		rowCount++;
	}

	delete[] temp_matrix;

	return 0;

}

int sparse::rowScale( sparse_matrix* A , int row1 , int row2 , double a ){

	/* create a new sparse matrix */
	sparse::sparse_matrix* temp_matrix = new sparse::sparse_matrix;
	temp_matrix->value_array = A->value_array;
	temp_matrix->colIdx = A->colIdx;
	temp_matrix->rowPtr = A->rowPtr;
	temp_matrix->noofRows = A->noofRows;
	temp_matrix->noofCols = A->noofCols;
	temp_matrix->noofVars = A->noofVars;

	/* Define some counters and limits */
	int nnCount_1 = A->rowPtr[row1];
	int nnLimit_1 = A->rowPtr[row1+1];
	int nnCount_2 = A->rowPtr[row2];
	int nnLimit_2 = A->rowPtr[row2+1];

	/* Find number of unique positions */
	int p = 0;
	bool* unique = new bool[temp_matrix->noofCols];

	for ( int i = 0 ; i < temp_matrix->noofCols ; i++ ) {
		unique[i] = false;
	}
	for ( int i = nnCount_1 ; i < nnLimit_1 ; i++ ) {
		p = temp_matrix->colIdx[i];
		unique[p] = true;
	}
	for ( int i = nnCount_2 ; i < nnLimit_2 ; i++ ) {
		p = temp_matrix->colIdx[i];
		unique[p] = true;
	}

	int unique_count = 0;
	for ( int i = 0 ; i < temp_matrix->noofCols ; i++ ) {
		if (unique[i]) unique_count++;
	}

	/* create new matrix */
	int sv = temp_matrix->noofVars - (nnLimit_2 - nnCount_2) + unique_count;
	A->value_array = new double[sv];
	A->colIdx = new int[sv];
	A->rowPtr = new int[temp_matrix->noofRows+1];
	A->noofVars = sv;

	/* copy elements upto row 2 */
	int count = 0;
	int rowCount = 0;
	for ( int i = 0 ; i < nnCount_2 ; i++ ) {
		A->value_array[count] = temp_matrix->value_array[i];
		A->colIdx[count] = temp_matrix->colIdx[i];
		count++;
	}

	for ( int i = 0 ; i <=row2 ; i++ ) {
		A->rowPtr[rowCount] = temp_matrix->rowPtr[i];
		rowCount++;
	}

	/* perform scaling operation on row 2 */
	for ( int i = 0 ; i < temp_matrix->noofCols ; i++ ) {
		if( unique[i] ){
			double temp , element1 , element2;
			sparse::retrieveElement(element1 , temp_matrix , row1 , i);
			sparse::retrieveElement(element2 , temp_matrix , row2 , i);
			temp = element2 + a*element1;
			A->value_array[count] = temp;
			A->colIdx[count] = i;
			count++;
		}
	}

	A->rowPtr[rowCount] = A->rowPtr[rowCount-1] + unique_count;
	rowCount++;

	/* copy elements from row 2 to end */
	for ( int i = nnLimit_2 ; i < temp_matrix->rowPtr[A->noofRows] ; i++ ){
		A->value_array[count] = temp_matrix->value_array[i];
		A->colIdx[count] = temp_matrix->colIdx[i];
		count++;
	}

	for ( int i = row2+1 ; i < temp_matrix->noofRows ; i++ ) {
		A->rowPtr[rowCount] = A->rowPtr[rowCount-1]+(temp_matrix->rowPtr[i+1] - temp_matrix->rowPtr[i]);
		rowCount++;
	}

	delete[] temp_matrix;
	delete[] unique;

	return 0;
}

int sparse::matrixProduct( double* result , sparse_matrix* A , double* vec ) {
	for ( int i = 0 ; i < A->noofRows ; i++ ){
		for ( int j = 0 ; j < A->noofCols ; j++ ) {
			double temp;
			retrieveElement( temp , A , i , j );
			result[i] += temp*vec[j];
		}
	}
	return 0;
}

int sparse::copyMatrix( sparse_matrix* CS , sparse_matrix* AS ) {
	/* create CS in main and pass pointer into this function */

	// copy over basic information, and allocate memory in heap
	CS->noofRows = AS->noofRows;
	CS->noofCols = AS->noofCols;
	CS->noofVars = AS->noofVars;
	CS->value_array = new double[CS->noofVars];
	CS->colIdx = new int[CS->noofVars];
	CS->rowPtr = new int[CS->noofRows+1];

	// copy over data
	for ( int i = 0 ; i < CS->noofVars ; i++ ){
		CS->value_array[i] = AS->value_array[i];
		CS->colIdx[i] = AS->colIdx[i];
	}

	for ( int i = 0 ; i < CS->noofRows+1 ; i++ ){
		CS->rowPtr[i] = AS->rowPtr[i];
	}

	return 0;
}

int sparse::decomposeMatrix( sparse_diagonal_matrix* DS , sparse_matrix* LUS , sparse_matrix* AS ) {
	
	/* create DS and LUS in main and pass pointer into this function */

	// create a copy of AS in LUS
	// extract the diagonal elements to DS
	// delete the diagonal elements of LUS

	copyMatrix( LUS , AS );
	DS->noofVars = AS->noofRows;
	DS->value_array = new double[DS->noofVars];

	for ( int i = 0 ; i < DS->noofVars ; i++ ){
		double temp;
		retrieveElement( temp , AS , i , i );
		DS->value_array[i] = temp;
		changeElement( LUS , i , i , 0.0 );
	}

	for ( int i = 0 ; i < LUS->noofRows ; i++ ){
		for ( int j = 0 ; j < LUS->noofCols ; j++ ){
			double temp;
			retrieveElement( temp , LUS , i , j );
			changeElement( LUS , i , j , -1.0*temp );
		}
	}

	return 0;
}

int sparse::jacobiSolver( double* X , sparse_diagonal_matrix* DS , sparse_matrix* LUS , double* B ){
	int size = LUS->noofRows;
	double* matPdt = new double[size];
	for ( int i = 0 ; i < size ; i++ ){
		matPdt[i] = 0.0;
	}
	matrixProduct( matPdt , LUS , X );

	for ( int i = 0 ; i < size ; i++ ){
		double dInv = 1.0/DS->value_array[i];
		X[i] = (B[i] + matPdt[i])*dInv;
	}
	return 0;
}