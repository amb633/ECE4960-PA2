#include <iostream>
#include <cstdlib>
#include <fstream>

using namespace std;

namespace sparse {
	typedef struct sparse_matrix {
		double* value_array;
		int* rowPtr;
		int* colIdx;
		int noofRows;
		int noofCols;
		int noofVars;
	} sparse_matrix;

	int retrieveElement( double& element , sparse_matrix* A , int rowInd , int colInd );
	int rowPermute ( sparse_matrix* A , int row1 , int row2 );
	int rowScale ( sparse_matrix* A , int row1 , int row2 , double a );
	int productAx( double* result , sparse_matrix* A , double* vec );

	int createMatrix( sparse_matrix* A );

}

namespace full{
	typedef struct full_matrix {
		double** value_array;
		double noofRows;
		double noofCols;
	} full_matrix;

	int retrieveElement( double& element , full_matrix* AF , int rowInd , int colInd );
	int rowPermute( full_matrix* AF , int row1 , int row2 );
	int rowScale( full_matrix* AF , int row1 , int row2 , double a );
	int productAx( double* result , full_matrix* AF , double* vec );
}

int main(int argc, char const *argv[])
{
	cout << endl;
	cout << "  /* Homework 04 : Programmed Individually by Haritha Murali (hm535) */ " << endl << endl;

	double* values = new double[12];
	for (int i = 0 ; i < 12 ; i++){
		values[i] = (double)(i+1);
	}

	int* row = new int[6];
	row[0] = 0;
	row[1] = 3;
	row[2] = 6;
	row[3] = 9;
	row[4] = 10;
	row[5] = 12;

	int* col = new int[12];
	col[0] = 0;
	col[1] = 1;
	col[2] = 4;
	col[3] = 0;
	col[4] = 1;
	col[5] = 2;
	col[6] = 1;
	col[7] = 2;
	col[8] = 4;
	col[9] = 3;
	col[10] = 0;
	col[11] = 4;

	sparse::sparse_matrix* A = new sparse::sparse_matrix;
	A->value_array = values;
	A->rowPtr = row;
	A->colIdx = col;
	A->noofRows = 5;
	A->noofCols = 5;
	A->noofVars = 12;

	full::full_matrix* AF = new full::full_matrix;
	AF->noofRows = 5;
	AF->noofCols = 5;

	AF->value_array = new double*[5];
	for ( int i = 0 ; i < 5 ; i++ ) {
		AF->value_array[i] = new double[5];
	}

	for ( int i = 0 ; i < 5 ; i++ ){
		for ( int j = 0 ; j < 5 ; j++ ){
			AF->value_array[i][j] = 0;
		}
	}

	AF->value_array[0][0] = 1;
	AF->value_array[0][1] = 2;
	AF->value_array[0][4] = 3;
	AF->value_array[1][0] = 4;
	AF->value_array[1][1] = 5;
	AF->value_array[1][2] = 6;
	AF->value_array[2][1] = 7;
	AF->value_array[2][2] = 8;
	AF->value_array[2][4] = 9;
	AF->value_array[3][3] = 10;
	AF->value_array[4][0] = 11;
	AF->value_array[4][4] = 12;

	cout << " Original Sparse Matrix : " << endl;

	for ( int i = 0 ; i < 5 ; i++ ) {
		cout << "   ";
		for ( int j = 0 ; j < 5 ; j++ ) {
			double element;
			sparse::retrieveElement( element , A , i , j );
			cout << element << "   " ;
		}
		cout << endl;
	}

	cout << endl;

	cout << " Original Full Matrix : " << endl;
	for ( int i = 0 ; i < 5 ; i++ ) {
		cout << "   ";
		for ( int j = 0 ; j < 5 ; j++ ) {
			double element;
			full::retrieveElement( element , AF , i , j );
			cout << element << "   ";
		}
		cout << endl;
	}

	cout << endl;

	cout << " Sparse Vector Product : " << endl;

	double v[5] = { 5 , 4 , 3 , 2 , 1 };
	double result_sparse[5] = { 0 , 0 , 0 , 0 , 0 };
	sparse::productAx( result_sparse , A , v );
	cout << "   " ;
	for ( int i = 0 ; i < 5 ; i++ ) {
		cout << result_sparse[i] << "   ";
	}
	cout << endl << endl;

	cout << " Full Vector Product : " << endl;
	double result_full[5] = { 0 , 0 , 0 , 0 , 0 };
	full::productAx( result_full , AF , v );
	cout << "   ";
	for ( int i = 0 ; i < 5 ; i++ ){
		cout << result_full[i] << "   ";
	}
	cout << endl << endl;

	cout << " 1st: Permute Rows 1 and 3 in Sparse Matrix : " << endl;

	sparse::rowPermute( A , 0 , 2 );

	for ( int i = 0 ; i < 5 ; i++ ) {
		cout << "   ";
		for ( int j = 0 ; j < 5 ; j++ ) {
			double element;
			sparse::retrieveElement( element , A , i , j );
			cout << element << "   " ;
		}
		cout << endl;
	}

	cout << endl << endl;

	cout << " 1st: Permute Rows 1 and 3 in Full Matrix : " << endl;

	full::rowPermute( AF , 0 , 2 );
	for ( int i = 0 ; i < 5 ; i++ ){
		cout << "   ";
		for ( int j = 0 ; j < 5 ; j++ ){
			double element;
			full::retrieveElement( element , AF , i , j );
			cout << element << "   ";
		}
		cout << endl;
	}

	cout << endl << endl;

	cout << " 2nd: Permute Rows 1 and 5 in Sparse Matrix : " << endl;

	sparse::rowPermute( A , 0 , 4 );

	for ( int i = 0 ; i < 5 ; i++ ) {
		cout << "   ";
		for ( int j = 0 ; j < 5 ; j++ ) {
			double element;
			sparse::retrieveElement( element , A , i , j );
			cout << element << "   " ;
		}
		cout << endl;
	}

	cout << endl << endl;

	cout << " 2nd: Permute Rows 1 and 5 in Full Matrix : " << endl;

	full::rowPermute( AF , 0 , 4 );
	for ( int i = 0 ; i < 5 ; i++ ){
		cout << "   ";
		for ( int j = 0 ; j < 5 ; j++ ){
			double element;
			full::retrieveElement( element , AF , i , j );
			cout << element << "   ";
		}
		cout << endl;
	}

	cout << endl << endl;


	cout << " Row Scale Operations ( with fresh matrices ) " << endl;

	/* First, create a new sparse matrix */
	sparse::sparse_matrix* B = new sparse::sparse_matrix;
	B->value_array = values;
	B->rowPtr = row;
	B->colIdx = col;
	B->noofRows = 5;
	B->noofCols = 5;
	B->noofVars = 12;

	/* Create new full matrix */
	full::full_matrix* BF = new full::full_matrix;
	BF->noofRows = 5;
	BF->noofCols = 5;

	BF ->value_array = new double*[5];
	for ( int i = 0 ; i < 5 ; i++ ){
		BF->value_array[i] = new double[5];
	}

	for ( int i = 0 ; i < 5 ; i++ ){
		for ( int j = 0 ; j < 5 ; j++ ){
			BF->value_array[i][j] = 0;
		}
	}

	BF->value_array[0][0] = 1;
	BF->value_array[0][1] = 2;
	BF->value_array[0][4] = 3;
	BF->value_array[1][0] = 4;
	BF->value_array[1][1] = 5;
	BF->value_array[1][2] = 6;
	BF->value_array[2][1] = 7;
	BF->value_array[2][2] = 8;
	BF->value_array[2][4] = 9;
	BF->value_array[3][3] = 10;
	BF->value_array[4][0] = 11;
	BF->value_array[4][4] = 12;

	cout << " 1st: Scale Rows 3.0*1 and 4 in Sparse Matrix : " << endl;

	sparse::rowScale( B , 0 , 3 , 3.0);

	for ( int i = 0 ; i < 5 ; i++ ) {
		cout << "   ";
		for ( int j = 0 ; j < 5 ; j++ ) {
			double element;
			sparse::retrieveElement( element , B , i , j );
			cout << element << "   " ;
		}
		cout << endl;
	}

	cout << endl << endl;

	cout << " 1sr: Scale Rows 3.0*1 and 4 in Full Matrix : " << endl;

	full::rowScale( BF , 0 , 3 , 3.0 );

	for ( int i = 0 ; i < 5 ; i++ ){
		cout << "   ";
		for ( int j = 0 ; j < 5 ; j++ ) {
			double element;
			full::retrieveElement( element , BF , i , j );
			cout << element << "   ";
		}
		cout << endl;
	}

	cout << endl << endl;

	cout << " 2nd: Scale Rows -4.4*5 and 2 in Sparse Matrix : " << endl;

	sparse::rowScale( B , 4 , 1 , -4.4);

	for ( int i = 0 ; i < 5 ; i++ ) {
		cout << "   ";
		for ( int j = 0 ; j < 5 ; j++ ) {
			double element;
			sparse::retrieveElement( element , B , i , j );
			cout << element << "   " ;
		}
		cout << endl;
	}

	cout << endl << endl;

	cout << " 2nd: Scale Rows -4.4*5 and 2 in Full Matrix : " << endl;

	full::rowScale( BF , 4 , 1 , -4.4 );

	for ( int i = 0 ; i < 5 ; i++ ){
		cout << "   ";
		for ( int j = 0 ; j < 5 ; j++ ){
			double element;
			full::retrieveElement( element , BF , i , j );
			cout << element << "   ";
		}
		cout << endl;
	}

	cout << endl << endl;

	cout << " ***** Creating First Large Matrix *****" << endl << endl;

	sparse::sparse_matrix* L1 = new sparse::sparse_matrix;
	sparse::createMatrix( L1 );

	cout << " Testing creation of Large Matrix : " << endl;
	double element;
	sparse::retrieveElement( element , L1 , 0 , 0 );
	cout << element << endl;
	sparse::retrieveElement( element , L1 , 1 , 0 );
	cout << element << endl;
	sparse::retrieveElement( element , L1 , 2 , 0 );
	cout << element << endl;
	sparse::retrieveElement( element , L1 , 2850 , 0 );
	cout << element << endl;

	int err = -1;
	err = sparse::rowPermute( L1 , 1 , 3 );
	if ( err == 0 ) cout << " first permute okay! " << endl;
	err = -1;
	err = sparse::rowPermute( L1 , 1 , 5 );
	if ( err == 0 ) cout << " second permute okay! " << endl;
	err = -1;
	err = sparse::rowPermute( L1 , 10 , 3000 );
	if ( err == 0 ) cout << " third permute okay! " << endl;
	err = -1;
	err = sparse::rowPermute( L1 , 5000 , 10000 );
	if ( err == 0 ) cout << " fourth permute okay! " << endl;
	err = -1;
	err = sparse::rowPermute( L1 , 6 , 15000 );
	if ( err == 0 ) cout << " fifth permute okay! " << endl;

	cout << endl<< " ***** Creating Second Large Matrix *****" << endl;

	sparse::sparse_matrix* L2 = new sparse::sparse_matrix;
	sparse::createMatrix( L2 );

	err = -1;
	err = sparse::rowScale( L2 , 2 , 4 , 3.0 );
	if ( err == 0 ) cout << " first scale okay! " << endl;
	err = -1;
	err = sparse::rowPermute( L2 , 2 , 5 );
	if ( err == 0 ) cout << " permute okay! " << endl;
	err = -1;
	err = sparse::rowScale( L2 , 5 , 4 , -3.0 );
	if ( err == 0 ) cout << " second scale okay! " << endl;

	cout << endl<< " ***** Creating Third Large Matrix *****" << endl;

	sparse::sparse_matrix* L3 = new sparse::sparse_matrix;
	sparse::createMatrix( L3 );

	/* Create vector */
	int size = L3->noofCols;
	double* vector_2 = new double[size];
	for ( int i = 0 ; i < size ; i++ ) {
		vector_2[i] = 1.0;
	}
	double* result_2 = new double[size];
	for ( int i = 0 ; i < size ; i++ ) {
		result_2[i] = 0.0;
	}
	err = -1;
	err = sparse::productAx( result_2 , L3 , vector_2 );
	if ( err == 0 ) cout << " matrix multiplication okay! " << endl;

	cout << endl << " ***** Calculating Difference in Sum of Elements ***** " << endl;
	/* Calculate the sum of elements in Large Matrix */

	cout.precision(10);

	double sum_matrix = 0.0;
	for ( int i = 0 ; i < L3->noofVars ; i++ ) {
		sum_matrix += L3->value_array[i];
	}

	cout << " sum of matrix elements = " << sum_matrix << endl;
	double sum_vector = 0.0;
	for ( int i = 0 ; i < L3->noofCols ; i++ ){
		sum_vector += result_2[i];
	}
	cout << " sum of vector product elements = " << sum_vector << endl;

	double error = sum_vector - sum_matrix;
	cout << " error between two sums = " << error << endl;

	cout << endl;
	return 0;
}

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

int sparse::productAx( double* result , sparse_matrix* A , double* vec ) {
	for ( int i = 0 ; i < A->noofRows ; i++ ){
		for ( int j = 0 ; j < A->noofCols ; j++ ) {
			double temp;
			retrieveElement( temp , A , i , j );
			result[i] += temp*vec[j];
		}
	}
	return 0;
}

int sparse::createMatrix( sparse_matrix* A ) {
	ifstream fin( "memplus.mtx" );
	int M , N , L;
	while ( fin.peek() == '%' ) fin.ignore( 2048 , '\n' );
	fin >> M >> N >> L;
	A->noofRows = M;
	A->noofCols = N;
	A->noofVars = L;
	A->value_array = new double[L];
	A->colIdx = new int[L];
	A->rowPtr = new int[M+1];

	int* keyed = new int[M+1];
	for ( int i = 0 ; i < M+1 ; i ++ ) {
		keyed[i] = 0;
	}

	for ( int i = 0 ; i < L ; i++ ) {
		int p;
		fin >> p;
		keyed[p]++;
		fin.ignore(2048, '\n');
	}

	for ( int i = 1 ; i < M+1 ; i++ ){
		keyed[i] += keyed[i-1];
	}

	fin.close();
	for ( int i = 0 ; i < M + 1 ; i++ ) {
		A->rowPtr[i] = keyed[i];
	}

	ifstream fil( "memplus.mtx" );
	while ( fil.peek() == '%' ) fil.ignore( 2048 , '\n' );
	fil.ignore( 2048 , '\n' );

	for ( int i = 0 ; i < L ; i++ ){
		int r , c;
		double data;
		fil >> r >> c >> data;
		int nnCount = keyed[r-1];
		A->value_array[nnCount] = data;
		A->colIdx[nnCount] = c-1;
		keyed[r-1]++;
	}

	fil.close();
	delete[] keyed;

	return 0;
}

int full::retrieveElement( double& element , full_matrix* AF , int rowInd , int colInd ) {
	element = AF->value_array[rowInd][colInd];
	return 0;
}

int full::productAx( double* result , full_matrix* AF , double* vec ){
	for ( int i = 0 ; i < 5 ; i++ ) {
		for ( int j = 0 ; j < 5 ; j++ ) {
			result[i] += AF->value_array[i][j]*vec[j];
		}
	}
	return 0;
}

int full::rowPermute( full_matrix* AF , int row1 , int row2 ){
	double* temp = AF->value_array[row2];
	AF->value_array[row2] = AF->value_array[row1];
	AF->value_array[row1] = temp;
	return 0;
}

int full::rowScale( full_matrix* AF , int row1 , int row2 , double a ){
	double* temp = AF->value_array[row1];
	for ( int i = 0 ; i < AF->noofCols ; i++ ){
		AF->value_array[row2][i] += temp[i]*a;
	}
	return 0;
}