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
   	
   	vector< vector<double> > input_vector = {{-4.0 , 1.0 , 0.0 , 0.0 , 1.0 }, {4.0 , -4.0 , 1.0 , 0.0 , 0.0 }, {0.0 , 1.0 , -4.0 , 1.0 , 0.0 }, {0, 0, 1, -4 , 1}, {1 , 0, 0, 1 , -4}};
   	vector< vector<double> > single_vector = {{1.0} , {0.0}, {0.0}, {0.0} , {0.0}};
   	
   	compressed::comp_r_mat mat_a = compressed::construct_compressed_matrix( &input_vector );

   	// const variable for for loops
   	const int RANK = mat_a.row_p.size() - 1;
   	cout << endl << "Original A matrix :" << endl;
   	for ( int i = 0 ; i < RANK ; i++ ){
    	for ( int j = 0 ; j < RANK ; j++ ){
           cout << compressed::retrieveElement( &mat_a , i , j ) << "   ";
       	}
       	cout << endl;
   	}

   	// testing copy function
   	cout << endl << "Copied C matrix : " << endl;
   	compressed::comp_r_mat mat_c;
   	compressed::copyMatrix( &mat_c , &mat_a );

   	for ( int i = 0 ; i < RANK ; i++ ){
    	for ( int j = 0 ; j < RANK ; j++ ){
        	cout << compressed::retrieveElement( &mat_c , i , j ) << "   ";
       	}
       	cout << endl;
   	}

   	// testing change non-zero element function
   	cout << endl << "C matrix with changed element : " << endl;
   	compressed::changeElement( &mat_c , 1 , 0 , 1.0 );

   	for ( int i = 0 ; i < RANK ; i++ ){
    	for ( int j = 0 ; j < RANK ; j++ ){
        	cout << compressed::retrieveElement( &mat_c , i , j ) << "   ";
       	}
       	cout << endl;
   	}

   	// testing if the copy made was indeed deep
   	cout << endl << "checking nothing was changed in A matrix (checking deep copy) " << endl;
   	for ( int i = 0 ; i < RANK ; i++ ){
       	for ( int j = 0 ; j < RANK ; j++ ){
           	cout << compressed::retrieveElement( &mat_a , i , j ) << "   ";
       	}
       	cout << endl;
   	}

   	// change back element in C matrix
   	compressed::changeElement( &mat_c , 1 , 0 , 4.0 );

   	// testing scalar multiplication
   	cout << endl << "checking scalar multiplication : " << endl;
   	compressed::scalarMultiple( &mat_c , 2.0 );
   	for ( int i = 0 ; i < RANK ; i++ ){
       	for ( int j = 0 ; j < RANK ; j++ ){
           	cout << compressed::retrieveElement( &mat_c , i , j ) << "   ";
       	}
       	cout << endl;
   	}

   	cout << endl;

   	// change back C matrix
   	compressed::scalarMultiple( &mat_c , 0.5 );
   
   	// testing matrix decomposition
   	compressed::comp_r_mat LUC;
   	compressed::comp_diag DC;
   	compressed::decomposeMatrix( &DC , &LUC , &mat_a );

   	cout << endl << "Diagonal Elements : " << endl;
   	for ( int i = 0 ; i < RANK ; i++ ){
    	cout << DC.value[i] << "   ";
   	}
   	cout << endl;

   	cout << endl << "LU Matrix : " <<  endl;
   	for ( int i = 0 ; i < RANK ; i++ ){
       	for ( int j = 0 ; j < RANK ; j++ ){
           	cout << compressed::retrieveElement( &LUC , i , j ) << "   ";
       	}
       	cout << endl;
   	}

   	compressed::comp_diag B; 
   	B.value.push_back(1.0); 
   	for ( int i = 0 ; i < RANK ; i++ ){
       	B.value.push_back(0.0);
   	}

   	cout << endl << "B Vector Elements : " << endl;
   	for ( int i = 0 ; i < RANK ; i++ ){
       	cout << B.value[i] << "   ";
   	}
   	cout << endl;

   	compressed::comp_diag X;
   	X.value.push_back(-0.25);
   	for (int i = 0; i < RANK ; i++ ) {
       	X.value.push_back(0.0);
   	}

   	cout << endl << "X Vector Elements : " << endl;
   	for ( int i = 0 ; i < RANK ; i++ ){
       	cout << X.value[i] << "   ";
   	}
   	cout << endl;

   	// testing matrix product
   	compressed::comp_diag mp;
   	for ( int i = 0 ; i < RANK ; i++ ){
       	mp.value.push_back(0.0);
   	}
   	compressed::matrixProduct( &mp , &mat_a , &B );
   	cout << endl << "mp Vector Elements : " << endl;
   	for ( int i = 0 ; i < RANK ; i++ ){
       	cout << mp.value[i] << "   ";
   	}
   	cout << endl << endl;

   	// reset mp to all zeros
   	for ( int i = 0 ; i < RANK ; i++ ){
   		mp.value[i] = 0.0;
   	}

   	// Jacobi Solver
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