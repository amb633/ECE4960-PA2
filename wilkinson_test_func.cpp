#include "wilkinson_test_func.hpp"

using namespace wilkinson_test;

bool wilkinson_test::test_retrieve_element( compressed::comp_r_mat* AC , vector< vector<double>>* AF ) {
	bool test = false; 
	int rank = (*AF).size();
	int counter = 0;
	for ( int i = 0 ; i < rank ; i++ ){
		for ( int j = 0 ; j < rank ; j++ ){
			if ( full::retrieveElement( AF,i,j ) == compressed::retrieveElement( AC,i,j ) )
				counter++;
		}
	}
	if ( counter == rank*rank ) test = true;
	return test;
}

bool wilkinson_test::test_change_element( compressed::comp_r_mat* AC , vector< vector<double>>* AF ) {
	bool test = false;
	double temp_c = compressed::retrieveElement( AC , 1 , 0 );
	double temp_f = full::retrieveElement( AF , 1 , 0 );
	compressed::changeElement( AC , 1 , 0 , 1.0 );
	full::changeElement( AF , 1 , 0 , 1.0 );
	if ( compressed::retrieveElement( AC , 1 , 0 ) == full::retrieveElement( AF , 1 , 0 ) )
		test = true;
	compressed::changeElement( AC , 1 , 0 , temp_c );
	full::changeElement( AF , 1 , 0 , temp_f );
	if ( compressed::retrieveElement( AC , 1 , 0 ) != full::retrieveElement( AF , 1 , 0 ) )
		test = false;
	return test;
}

bool wilkinson_test::test_copy_matrix( compressed::comp_r_mat* AC , vector< vector<double>>* AF ){
	bool test = false;
	compressed::comp_r_mat AC_copy;
	vector< vector<double>> AF_copy;
	compressed::copyMatrix( &AC_copy , AC );
	full::copyMatrix( &AF_copy , AF );
	if (test_retrieve_element( &AC_copy , &AF_copy ) ) test = true;
	return test;
}

bool wilkinson_test::test_scalar_multiple( compressed::comp_r_mat* AC , vector< vector<double>>* AF ){
	bool test = false;
	compressed::scalarMultiple( AC , 2.0 );
	full::scalarMultiple( AF , 2.0 );
	if ( test_retrieve_element ( AC , AF ) ) test = true;
	compressed::scalarMultiple( AC , 0.5 );
	full::scalarMultiple( AF , 0.5 );
	if ( !test_retrieve_element(AC , AF) ) test = false;
	return test;
}

bool wilkinson_test::test_row_permute( compressed::comp_r_mat* AC , vector< vector<double>>* AF ){
	bool test = false;
	compressed::rowPermute( AC , 1 , 3 );
	full::rowPermute( AF , 1 , 3 );
	if ( test_retrieve_element ( AC , AF ) ) test = true;
	compressed::rowPermute( AC , 1 , 3 );
	full::rowPermute( AF , 1 , 3 );
	if ( !test_retrieve_element( AC , AF ) ) test = false;
	return test;
}

bool wilkinson_test::test_row_scale( compressed::comp_r_mat* AC , vector< vector<double>>* AF ){
	bool test = false;
	compressed::rowScale( AC , 0 , 2 , 3.0 );
	full::rowScale( AF , 0 , 2 , 3.0 );
	if ( test_retrieve_element ( AC , AF ) ) test = true;
	compressed::rowScale( AC , 0 , 2 , -3.0 );
	full::rowScale( AF , 0 , 2 , -3.0 );
	if ( !test_retrieve_element( AC , AF ) ) test = false;
	return test;
}

bool wilkinson_test::test_matrix_product( compressed::comp_r_mat* AC , vector< vector<double>>* AF , vector<double>* B ){
	bool test = false;
	int rank = (*B).size();
	vector<double> result_comp , result_full;
	for ( int i = 0 ; i < rank ; i++ ){
		result_comp.push_back(0.0);
		result_full.push_back(0.0);
	}
	compressed::productAx( &result_comp , AC , B );
	full::productAx( &result_full , AF , B );

	int counter = 0;
	for ( int i = 0 ; i < rank ; i++ ){
		if ( result_comp[i] == result_full[i] ) counter++;
	}

	if ( counter == rank ) test = true;
	return test;

}

bool wilkinson_test::test_calculate_norm( compressed::comp_r_mat* AC , vector< vector<double>>* AF , vector<double>* B ){
	bool test = false;
	int rank = (*B).size();
	vector<double> result_comp , result_full;
	for ( int i = 0 ; i < rank ; i++ ){
		result_comp.push_back(0.0);
		result_full.push_back(0.0);
	}
	compressed::productAx( &result_comp , AC , B );
	full::productAx( &result_full , AF , B );

	double norm_comp , norm_full;
	compressed::calculateNorm( norm_comp , &result_comp , &result_full );
	full::calculateNorm( norm_full , &result_full , &result_comp );

	if ( norm_full == norm_comp ) test = true;
	return test;

}

bool wilkinson_test::test_matrix_decomposition( compressed::comp_r_mat* AC , vector< vector<double>>* AF ){
	bool test = false;
	int rank = (*AF).size();
	int counter = 0;

	vector<double> DC , DF;
	vector< vector<double>> LUF;
	compressed::comp_r_mat LUC;

	compressed::decomposeMatrix( &DC , &LUC , AC );
	full::decomposeMatrix( &DF , &LUF , AF );

	if ( test_retrieve_element ( &LUC , &LUF ) ) counter++;
	for ( int i = 0 ; i < rank ; i++ ){
		if ( DC[i] == DF[i] ) counter++;
	}
	if ( counter == rank + 1 ) test = true;
	return test;
}

bool wilkinson_test::test_jacobi_solver( compressed::comp_r_mat* AC , vector< vector<double>>* AF  , vector<double>* B ){
	bool test = false;
	int rank = (*B).size();

	vector<double> X_F , X_C;
	vector<double> DC , DF;
	vector<double> mpc , mpf;

	vector< vector<double>> LUF;
	compressed::comp_r_mat LUC;

   	X_F.push_back(-0.25);
   	X_C.push_back(-0.25);
   	mpc.push_back(0.0);
	mpf.push_back(0.0);

   	for (int i = 0; i < rank - 1 ; i++ ) {
       	X_F.push_back(0.0);
       	X_C.push_back(0.0);
       	mpc.push_back(0.0);
       	mpf.push_back(0.0);
   	}

   	compressed::decomposeMatrix( &DC , &LUC , AC );
   	double normPrev = 2;
    double normCurrent = 1;
    while ( abs(normCurrent - normPrev) > 1e-7) {
       	normPrev = normCurrent;
       	compressed::jacobiSolver( &X_C , &DC , &LUC , B );
       	compressed::productAx( &mpc , AC , &X_C );
       	compressed::calculateNorm( normCurrent , B , &mpc );
   	} 

   	full::decomposeMatrix( &DF , &LUF , AF );
   	normPrev = 2;
   	normCurrent = 1;
   	while ( abs(normCurrent - normPrev) > 1e-7) {
       	normPrev = normCurrent;
       	full::jacobiSolver( &X_F , &DF , &LUF , B );
       	full::productAx( &mpf , AF , &X_F );
       	full::calculateNorm( normCurrent , B , &mpf );
   	} 

   	int counter = 0;
   	for ( int i = 0 ; i < rank ; i++ ){
   		if ( X_C[i] == X_F[i] ) counter++;
   	}
   	if ( counter == rank ) test = true;
   	return test;
}

int main(int argc, char const *argv[])
{
	cout << endl << boolalpha;

	vector< vector<double>> A = {{-4.0 , 1.0 , 0.0 , 0.0 , 1.0 }, {4.0 , -4.0 , 1.0 , 0.0 , 0.0 }, {0.0 , 1.0 , -4.0 , 1.0 , 0.0 }, {0, 0, 1, -4 , 1}, {1 , 0, 0, 1 , -4}};
	vector<double> B = {1.0 , 0.0 , 0.0 , 0.0 , 0.0 };

	compressed::comp_r_mat AC = compressed::construct_compressed_matrix( &A );

	cout << "testing retrieve element functions : " << test_retrieve_element( &AC , &A ) << endl;
	cout << "testing change element functions   : " << test_change_element( &AC , &A ) << endl;
	cout << "testing copy matrix functions      : " << test_copy_matrix( &AC , &A ) << endl;
	cout << "testing scalar multiple functions  : " << test_scalar_multiple( &AC , &A ) << endl;
	cout << "testing row permutation functions  : " << test_row_permute( &AC , &A ) << endl;
	cout << "testing row scaling functions      : " << test_row_scale( &AC , &A ) << endl;
	cout << "testing matrix product functions   : " << test_matrix_product( &AC , &A , &B ) << endl;
	cout << "testing calculate norm functions   : " << test_calculate_norm( &AC , &A , &B ) << endl;
	cout << "testing matrix decomposition funcs : " << test_matrix_decomposition( &AC , &A ) << endl;
	
	cout << "testing jacobi solver functions #1 : " << test_jacobi_solver( &AC , &A , &B ) << endl;

	B[0] = 0.0;
	B[3] = 1.0;

	cout << "testing jacobi solver functions #2 : " << test_jacobi_solver( &AC , &A , &B ) << endl;

	for ( int i = 0 ; i < B.size() ; i++ ){
		B[i] = 1.0;
	}

	cout << "testing jacobi solver functions #3 : " << test_jacobi_solver( &AC , &A , &B ) << endl;

	compressed::rowScale( &AC , 0 , 2 , 3.0 );
	full::rowScale( &A , 0 , 2 , 3.0 );

	/*for ( int i = 0 ; i < B.size() ; i++ ){
		for ( int j = 0 ; j < B.size() ; j++ ){
			cout << compressed::retrieveElement( &AC , i , j ) << "   ";
		}
		cout << endl;
	}
	cout << endl << endl;
	for ( int i = 0 ; i < B.size() ; i++ ){
		for ( int j = 0 ; j < B.size() ; j++ ){
			cout << full::retrieveElement( &A , i , j ) << "   ";
		}
		cout << endl;
	}*/



	cout << endl;
	return 0;
}

