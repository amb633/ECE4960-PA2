#include "large_matrix_test.hpp"

using namespace large_matrix_test;

bool large_matrix_test::test_row_permute( compressed::comp_r_mat* AC ){
	bool test = true;
	int rank = (*AC).row_p.size() - 1;
	// 1: create a copy of the matrix
	compressed::comp_r_mat CC;
	compressed::copyMatrix( &CC , AC );
	//2: permute rows in the copied matrix
	compressed::rowPermute( &CC , 1500 , 2500 );
	//3: test if elements in row1 of copy = elements in row2 of original
	for ( int i = 0 ; i < rank ; i++ ){
		if ( compressed::retrieveElement( &CC , 1500 , i ) != compressed::retrieveElement( AC , 2500 , i ) ) 
			test = false;
	}
	//4: rest if elements in row2 to copy = elements in row1 of original
	for (int i = 0 ; i < rank ; i++ ){
		if ( compressed::retrieveElement( &CC , 2500 , i ) != compressed::retrieveElement( AC , 1500 , i ) )
			test = false;
	}
	return test;
}

bool large_matrix_test::test_row_scale( compressed::comp_r_mat* AC ){
	bool test = true;
	int rank = (*AC).row_p.size() - 1;
	//1: create a copy of the matrix
	compressed::comp_r_mat CC;
	compressed::copyMatrix( &CC , AC );
	//2: scale rows in copied matrix
	compressed::rowScale( &CC , 500 , 2000 , 2.0 );
	//3: test elements in row 2 of copy
	for( int i = 0 ; i < rank ; i++ ){
		if ( compressed::retrieveElement( &CC , 2000 , i ) != 2.0*compressed::retrieveElement( AC , 500 , i ) + compressed::retrieveElement( AC , 2000 , i ) )
			test = false;
	}
	//4: test elements in row 1 of copy
	for( int i = 0 ; i < rank ; i++ ){
		if ( compressed::retrieveElement( &CC , 500 , i ) != compressed::retrieveElement( AC , 500 , i ) )
			test = false;
	}
	return test;
}

bool large_matrix_test::test_matrix_product( compressed::comp_r_mat* AC ){
	bool test = true;
	int rank = (*AC).row_p.size() - 1;
	//1: create vectors for multiplication
	vector<double> product , B;
	for ( int i = 0 ; i < rank ; i++ ){
		B.push_back(1.0);
		product.push_back(0.0);
	}
	//2: call the matrix product function
	compressed::productAx( &product , AC , &B );
	//3: check sum
	double sum_vector = 0.0 ;
	double sum_matrix = 0.0 ;
	for ( int i = 0 ; i < rank ; i++ ){
		sum_vector += product[i];
	}
	for ( int i = 0 ; i < (*AC).value.size() ; i++ ){
		sum_matrix += (*AC).value[i];
	}
	if ( abs(sum_matrix - sum_vector) > 10e-8 ) test = false;
	return test;
}

bool large_matrix_test::test_calculate_norm( compressed::comp_r_mat* AC ){
	bool test = false;
	int rank = (*AC).row_p.size() - 1;
	//1: create a copy of AC
	compressed::comp_r_mat CC;
	compressed::copyMatrix( &CC , AC );
	//2: permute some rows in CC
	compressed::rowPermute( &CC , 500 , 2000 );
	compressed::rowPermute( &CC , 1500 , 3000 );
	compressed::rowPermute( &CC , 4000 , 4500 );
	//3: create vectors for multiplication
	vector<double> product_A , product_C , B;
	for ( int i = 0 ; i < rank ; i++ ){
		product_A.push_back(0.0);
		product_C.push_back(0.0);
		B.push_back(1.0);
	}
	//4: call matrix product functions
	compressed::productAx( &product_A , AC , &B );
	compressed::productAx( &product_C , &CC , &B );
	//5: reverse-permute the rows in product_C
	double temp = product_C[500];
	product_C[500] = product_C[2000];
	product_C[2000] = temp;
	temp = product_C[1500];
	product_C[1500] = product_C[3000];
	product_C[3000] = temp;
	temp = product_C[4000];
	product_C[4000] = product_C[4500];
	product_C[4500] = temp;
	//6: the vectors should be identical
	double norm;
	compressed::calculateNorm( norm , &product_C , &product_A );
	if ( norm < 10e-11 ) test = true;
	return test;
}

bool large_matrix_test::test_matrix_decomposition( compressed::comp_r_mat* AC ){
	bool test = true;
	int rank = (*AC).row_p.size() - 1;
	//1: create variables for storing decomposed matrix
	vector<double> diagonals;
	compressed::comp_r_mat LUC;
	//2: call the decomposition function
	compressed::decomposeMatrix( &diagonals , &LUC , AC );
	//3: check if the diagonals are correct
	for ( int i = 0 ; i < rank ; i++ ){
		if ( diagonals[i] != compressed::retrieveElement( AC , i , i ) ) test = false;
		if( compressed::retrieveElement( &LUC , i , i ) != 0.0 ) test = false;
	}
	return test;
}

int main(int argc, char const *argv[])
{
	cout << endl;
	ifstream rowPtr_file("./Mat1/rowPtr.csv");
    int row_value;
    vector<int> row_ptr;
    while(rowPtr_file >> row_value){
        row_ptr.push_back(row_value-1);
    }
    cout << "This is the number of row_ptrs in Mat 1: " << row_ptr.size() << endl;
    
    ifstream values_file("./Mat1/value.csv");
    double value;
    vector<double> values;
    while(values_file >> value){
        values.push_back(value);
    }
    cout << "This is the number of non_zero values in Mat 1: " << values.size() << endl;
    
    ifstream col_file("./Mat1/colInd.csv");
    int col_id;
    vector<int> cols;
    while(col_file >> col_id){
        cols.push_back(col_id-1);
    }
    cout << "This is the number of col_ids for non-zero values in Mat 1: " << cols.size() << endl;

    const int rank = row_ptr.size() - 1;

    compressed::comp_r_mat AC;
    AC.value = values;
    AC.col_id = cols;
    AC.row_p = row_ptr;

    cout << endl << boolalpha;
    cout << "testing large matrix row permutation: " << test_row_permute( &AC ) << endl;
    //cout << "testing large matrix row scaling    : " << test_row_scale ( &AC ) << endl;
    cout << "testing large matrix product        : " << test_matrix_product( &AC ) << endl;
    cout << "testing large matrix calculate norm : " << test_calculate_norm( &AC ) << endl;
    cout << "testing large matrix decomposition  : " << test_matrix_decomposition( &AC ) << endl;

    cout << endl;
	return 0;
}