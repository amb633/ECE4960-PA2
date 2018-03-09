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
   	compressed::comp_diag DC;
   	compressed::decomposeMatrix( &DC , &LUC , &mat_a );

   	cout << endl << "Diagonal Elements : Rank = " << DC.rank << endl;
   	for ( int i = 0 ; i < DC.rank ; i++ ){
    	cout << DC.value[i] << "   ";
   	}
   	cout << endl;

   	cout << endl << "LU Matrix : Rank Info : " << LUC.noofRows << " , " << LUC.noofCols << " , " << LUC.noofVars << endl;
   	for ( int i = 0 ; i < 5 ; i++ ){
       	for ( int j = 0 ; j < 5 ; j++ ){
           	cout << compressed::retrieveElement( &LUC , i , j ) << "   ";
       	}
       	cout << endl;
   	}

   	compressed::comp_diag B;
   	B.rank = 5; 
   	B.value.push_back(1.0); 
   	for ( int i = 0 ; i < B.rank ; i++ ){
       	B.value.push_back(0.0);
   	}

   	cout << endl << "B Vector Elements : Rank = " << B.rank << endl;
   	for ( int i = 0 ; i < B.rank ; i++ ){
       	cout << B.value[i] << "   ";
   	}
   	cout << endl;

   	compressed::comp_diag X;
   	X.rank = 5;
   	X.value.push_back(-0.25);
   	for (int i = 0; i < X.rank ; i++ ) {
       	X.value.push_back(0.0);
   	}

   	cout << endl << "X Vector Elements : Rank = " << X.rank << endl;
   	for ( int i = 0 ; i < X.rank ; i++ ){
       	cout << X.value[i] << "   ";
   	}
   	cout << endl;

   	// testing matrix product
   	compressed::comp_diag mp;
   	mp.rank = B.rank;
   	for ( int i = 0 ; i < 5 ; i++ ){
       	mp.value.push_back(0.0);
   	}
   	compressed::matrixProduct( &mp , &mat_a , &B );
   	cout << endl << "mp Vector Elements : Rank = " << mp.rank << endl;
   	for ( int i = 0 ; i < mp.rank ; i++ ){
       	cout << mp.value[i] << "   ";
   	}
   	cout << endl;

   	// reset mp to all zeros
   	for ( int i = 0 ; i < mp.rank ; i++ ){
   		mp.value[i] = 0.0;
   	}

    double normPrev = 2;
    double normCurrent = 1;
    int counter = 0;

  	// for ( int j = 0 ; j < 50 ; j++ ) {
  	while ( abs(normCurrent - normPrev) > 1e-7) {
       	normPrev = normCurrent;
       	compressed::jacobiSolver( &X , &DC , &LUC , &B );
       	compressed::matrixProduct( &mp , &mat_a , &X );
       	compressed::calculateNorm( normCurrent , &B , &mp );

       	cout << counter << " : " ;
       	for ( int i = 0 ; i < 5 ; i++ ){
           	cout << X.value[i] << "   ";
       	}
       	cout << " : " << normCurrent << endl;
       	counter++;
   	}

   	cout<<endl;
   	return 0;
}