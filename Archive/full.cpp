#include "full.hpp"

int full::retrieveElement( double& element , full_matrix* AF , int rowInd , int colInd ) {
	element = AF->value_array[rowInd][colInd];
	return 0;
}

int full::matrixProduct( double* result , full_matrix* AF , double* vec ){
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

int full::copyMatrix( full_matrix* CF , full_matrix* AF ){
	/* create CF in main and pass pointer into this function */
	int size = AF->noofRows;
	CF->noofRows = size;
	CF->noofCols = size;

	CF->value_array = new double*[size];
	for ( int i = 0 ; i < size ; i++ ){
		CF->value_array[i] = new double[size];
	}

	for( int i = 0 ; i < CF->noofRows ; i++ ){
		for ( int j = 0 ; j < CF->noofCols ; j++ ){
			CF->value_array[i][j] = AF->value_array[i][j];
		}
	}
}

int full::decomposeMatrix( full_diagonal_matrix* DF , full_matrix* LUF , full_matrix* AF ){
	
	copyMatrix( LUF , AF );
	DF->noofVars = AF->noofRows;
	DF->value_array = new double[DF->noofVars];
	
	for ( int i = 0 ; i < DF->noofVars ; i++ ){
		double temp;
		retrieveElement( temp , AF , i , i );
		DF->value_array[i] = temp;
		LUF->value_array[i][i] = 0.0;
	}
}