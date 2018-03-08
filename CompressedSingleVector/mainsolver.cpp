#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>

#include "compressed_func.hpp"
#include "full_func.hpp"

using namespace std;

int main(int argc, char const *argv[])
{
	// testing creation of compressed matrix
	cout << endl << "Original A matrix :" << endl;
	vector< vector<double> > input_vector = {{-4 , 1, 0, 0 , 1}, {4 , -4, 1, 0, 0}, {0 , 1, -4, 1, 0}, {0, 0, 1, -4 , 1}, {1 , 0, 0, 1 , -4}};
	vector< vector<double> > single_vector = {{1} , {0}, {0}, {0} , {0}};
	
	compressed::comp_r_mat mat_a = compressed::construct_compressed_matrix( &input_vector );
	mat_a.noofRows = 5;
	mat_a.noofCols = 5;
	mat_a.noofVars = 15;

	for ( int i = 0 ; i < 5 ; i++ ){
		for ( int j = 0 ; j < 5 ; j++ ){
			cout << compressed::retrieveElement( &mat_a , i , j ) << "   ";
		}
		cout << endl;
	}

	// testing copy function
	cout << endl << "Copied C matrix : " << endl;
	compressed::comp_r_mat mat_c;
	compressed::copyMatrix( &mat_c , &mat_a );

	for ( int i = 0 ; i < 5 ; i++ ){
		for ( int j = 0 ; j < 5 ; j++ ){
			cout << compressed::retrieveElement( &mat_c , i , j ) << "   ";
		}
		cout << endl;
	}

	cout << " extra info: " << mat_c.noofRows << " , " << mat_c.noofCols << " , " << mat_c.noofVars << endl;

	// testing change nz element function
	cout << endl << "C matrix with changed element : " << endl;
	compressed::changeElement( &mat_c , 1 , 0 , 1.0 );

	for ( int i = 0 ; i < 5 ; i++ ){
		for ( int j = 0 ; j < 5 ; j++ ){
			cout << compressed::retrieveElement( &mat_c , i , j ) << "   ";
		}
		cout << endl;
	}

	// testing if the copy made was indeed deep
	cout << endl << "checking nothing was changed in A matrix (checking deep copy) " << endl;
	for ( int i = 0 ; i < 5 ; i++ ){
		for ( int j = 0 ; j < 5 ; j++ ){
			cout << compressed::retrieveElement( &mat_a , i , j ) << "   ";
		}
		cout << endl;
	}

	// change back element in C matrix
	compressed::changeElement( &mat_c , 1 , 0 , 4.0 );

	// testing scalar multiplication
	cout << endl << "checking scalar multiplication : " << endl;
	compressed::scalarMultiple( &mat_c , -1.0 );
	for ( int i = 0 ; i < 5 ; i++ ){
		for ( int j = 0 ; j < 5 ; j++ ){
			cout << compressed::retrieveElement( &mat_c , i , j ) << "   ";
		}
		cout << endl;
	}

	cout << endl;

	// change back C matrix
	compressed::scalarMultiple( &mat_c , -1.0 );
	
	// testing matrix decomposition
	compressed::comp_r_mat LUC;
	compressed::comp_r_mat DC;
	compressed::decomposeMatrix( &DC , &LUC , &mat_a );

	cout << endl << "Diagonal Elements : Rank = " << DC.noofRows << endl;
	for ( int i = 0 ; i < DC.noofRows ; i++ ){
		cout << compressed::retrieveElement(&DC , i , i ) << "   ";
	}
	cout << endl;

	cout << endl << "LU Matrix : Rank Info : " << LUC.noofRows << " , " << LUC.noofCols << " , " << LUC.noofVars << endl;
	for ( int i = 0 ; i < 5 ; i++ ){
		for ( int j = 0 ; j < 5 ; j++ ){
			cout << compressed::retrieveElement( &LUC , i , j ) << "   ";
		}
		cout << endl;
	}

	compressed::comp_r_mat mat_b2;
	mat_b2.noofRows = 5;
	mat_b2.noofCols = 1;
	mat_b2.noofVars = 1;

	int row = 0;
    for( int i = 0; i< mat_b2.noofRows; i++){
        mat_b2.row_p.push_back(row);
        if( i == 0 ) {
            mat_b2.value.push_back(1.0);
            mat_b2.col_id.push_back(0);
            row++;
        }
        mat_b2.row_p.push_back(row);
    }

    cout << endl << "B Vector Elements : Rank = " << mat_b2.noofRows << endl;
	for ( int i = 0 ; i < mat_b2.noofRows ; i++ ){
		cout << compressed::retrieveElement(&mat_b2 , i , 0 ) << "   ";
	}
	cout << endl;

	compressed::comp_r_mat X;
	X.noofRows = 5;
	X.noofCols = 1;
	X.noofVars = 5;

	row = 0;
	for ( int i = 0 ; i < X.noofRows ; i++ ){
		X.row_p.push_back(row);
		if ( i == 0 ){
			X.value.push_back(-0.25);
			X.col_id.push_back(0);
			row++;
		}
		X.row_p.push_back(row);
	}

	cout << endl << "X Vector Elements : Rank = " << X.noofRows << endl;
	for ( int i = 0 ; i < X.noofRows ; i++ ){
		cout << compressed::retrieveElement(&X , i , 0 ) << "   ";
	}
	cout << endl;

	// testing matrix product
	compressed::comp_r_mat product;
	product.noofRows = 5;
	product.noofCols = 1;
	product.noofVars = 5;

	compressed::matrixProduct( &product , &mat_a , &mat_b2 );
	cout << endl << "matrix product : " << endl;
	for ( int i = 0 ; i < product.noofRows ; i++ ){
		cout << compressed::retrieveElement( &product , i , 0 ) << "   ";
	}
	cout << endl << endl ;

	cout << endl << "jacobi iterations : " << endl;
	int counter = 0;
	for ( int j = 0 ; j < 5 ; j++ )
	{
		counter++;
		compressed::jacobiSolver( &X , &DC , &LUC , &mat_b2 );
		//cout << counter << " : " ;
		for ( int i = 0 ; i < 5 ; i++ ){
			//cout << compressed::retrieveElement( &X , i , 0 ) << "   ";
		}
		cout << endl;
	}

	cout<<endl;
	return 0;
}