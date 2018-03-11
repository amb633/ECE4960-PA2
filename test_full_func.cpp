#include "test_full_func.hpp"

vector< vector<double> > AF = { { 1 , 2 , 0 , 0 , 3 } , { 4 , 5 , 6 , 0 , 0 } , { 0 , 7 , 8 , 0 , 9 } , { 0 , 0 , 0 , 10 , 0 } , { 11 , 0 , 0 , 0 , 12 } };
vector<double> VF = { 5 , 4 , 3 , 2 , 1 };

bool test_full::test_matrix_retrieve_element(){
	bool test = false;
	double element;
	vector<double> ground_truth = { 1 , 5 , 8 , 10 , 12 };

	int counter = 0;
	for ( int i = 0 ; i < AF.size() ; i++ ){
		element = retrieveElement( &AF , i , i );
		if ( element == ground_truth[i] ) counter++;
	}
	if ( counter == ground_truth.size() ) test = true;

	return test;
}

bool test_full::test_vector_retrieve_element(){
	bool test = false;
	double element;
	vector<double> ground_truth = { 5 , 4 , 3 , 2 , 1 };

	int counter = 0;
	for ( int i = 0 ; i < VF.size() ; i++ ){
		element = retrieveElement( &VF , i );
		if ( element == ground_truth[i] ) counter++;
	}
	if ( counter == ground_truth.size() ) test = true;

	return test;
}

bool test_full::test_change_element(){
	bool test = false;
	changeElement( &AF , 3 , 1 , 15.0 );
	if ( AF[3][1] == 15.0 ) test = true;
	changeElement( &AF , 3 , 1 , 0.0 );
	if ( AF[3][1] != 0.0 ) test = false;
	return test;
}

bool test_full::test_copy_matrix(){
	bool test = false;
	vector< vector<double>> CF;
	double element;

	copyMatrix( &CF , &AF );

	// test if the copied matrix is identical
	int counter = 0;
	for ( int i = 0 ; i < AF.size() ; i++ ){
		for ( int j = 0 ; j < AF.size() ; j++ ){
			if ( CF[i][j] == AF[i][j] )	counter++;
		}
	}
	// test if indeed a deep copy
	changeElement( &CF , 3 , 1 , 15.0 );
	if ( AF[3][1] == 0.0 ) counter++;
	if ( counter == (AF.size()*AF.size() + 1) ) test = true;
	return test;
}

bool test_full::test_scalar_multiple(){
	bool test = false;
	vector< vector<double>> CF;
	copyMatrix( &CF , &AF );
	scalarMultiple( &CF , 2.0 );
	int counter = 0;
	for ( int i = 0 ; i < AF.size() ; i++ ){
		for ( int j = 0 ; j < AF.size() ; j++ ){
			if ( CF[i][j] == 2*AF[i][j] ) counter++;
		}
	}
	if ( counter == (AF.size()*AF.size()) ) test = true;
	return test;
}

bool test_full::test_row_permute(){
	bool test = false;

	vector< vector<double>> CF;
	copyMatrix( &CF , &AF );
	rowPermute( &CF , 0 , 3 );
	int counter = 0;
	for ( int i = 0 ; i < AF.size() ; i++ ){
		if ( CF[0][i] == AF[3][i] ) counter++;
		if ( CF[1][i] == AF[1][i] ) counter++;
		if ( CF[2][i] == AF[2][i] ) counter++;
		if ( CF[3][i] == AF[0][i] ) counter++;
		if ( CF[4][i] == AF[4][i] ) counter++;
	}
	if ( counter == (AF.size()*AF.size()) ) test = true;
	return test;

}

bool test_full::test_row_scale(){
	bool test = false;
	vector< vector<double>> CF;
	copyMatrix( &CF , &AF );
	rowScale( &CF , 0 , 4 , 3.0 );
	
	int counter = 0;
	for ( int i = 0 ; i < AF.size() - 1 ; i++ ){
		for ( int j = 0 ; j < AF.size() ; j++ ){
			if ( CF[i][j] == AF[i][j] ) counter++;
		}
	}
	for ( int i = 0 ; i < AF.size() ; i++ ){
		if( CF[4][i] == 3.0*AF[0][i] + AF[4][i] ) counter++;
	}

	if ( counter == (AF.size()*AF.size()) ) test = true;
	return test;
}

bool test_full::test_matrix_product(){
	bool test = false;
	vector<double> prod;
	for ( int i = 0 ; i < VF.size() ; i++ ){
		prod.push_back( 0.0 );
	}

	vector<double> ground_truth = { 16.0 , 58.0 , 61.0 , 20.0 , 67.0 };
	productAx( &prod , &AF , &VF );

	int counter = 0;
	for ( int i = 0 ; i < prod.size() ; i++ ){
		if( prod[i] == ground_truth[i] ) counter++;
	}
	if ( counter == VF.size() ) test = true;
	return test;
}

bool test_full::test_calculate_norm(){
	bool test = false;
	vector<double> prod;
	vector<double> zeros = { 0.0 , 0.0 , 0.0 , 0.0 , 0.0 };
	vector<double> ones = { 1.0 , 1.0 , 0.0 , -1.0 , -1.0 };
	vector<double> twos = { 2.0 , 2.0 , 0.0 , -2.0 , -2.0 };

	double norm_ones , norm_twos;
	calculateNorm( norm_ones , &zeros , &ones );
	calculateNorm( norm_twos , &zeros , &twos );
	if( norm_ones == 2.0 && norm_twos == 4.0 ) test = true;
	return test;
}

bool test_full::test_matrix_decomposition(){
	bool test = false;
	vector<double> diagonals = { 1 , 5 , 8 , 10 , 12 };
	vector< vector<double> > LU_true = { { 0 , -2 , 0 , 0 , -3 } , { -4 , 0 , -6 , 0 , 0 } , { 0 , -7 , 0 , 0 , -9 } , { 0 , 0 , 0 , 0 , 0 } , { -11 , 0 , 0 , 0 , 0 } };
	vector<double> DF;
	vector< vector<double>> LUF;
	decomposeMatrix( &DF , &LUF , &AF );
	int counter = 0;
	for ( int i = 0 ; i < AF.size() ; i++ ){
		for ( int j = 0 ; j < AF.size() ; j++ ){
			if ( LUF[i][j] == LU_true[i][j] ) counter++;
		}
		if( DF[i] == diagonals[i] ) counter++;
	}
	if ( counter == ( (AF.size()*AF.size()) + AF.size() ) ) test = true;
	return test;
}