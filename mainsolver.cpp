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
   	vector<double> DC;
   	compressed::decomposeMatrix( &DC , &LUC , &mat_a );

   	cout << endl << "Diagonal Elements : " << endl;
   	for ( int i = 0 ; i < RANK ; i++ ){
    	cout << DC[i] << "   ";
   	}
   	cout << endl;

   	cout << endl << "LU Matrix : " <<  endl;
   	for ( int i = 0 ; i < RANK ; i++ ){
       	for ( int j = 0 ; j < RANK ; j++ ){
           	cout << compressed::retrieveElement( &LUC , i , j ) << "   ";
       	}
       	cout << endl;
   	}

   	vector<double> B; 
   	B.push_back(1.0); 
   	for ( int i = 0 ; i < RANK - 1 ; i++ ){
       	B.push_back(0.0);
   	}

   	cout << endl << "B Vector Elements : " << endl;
   	for ( int i = 0 ; i < B.size() ; i++ ){
       	cout << B[i] << "   ";
   	}
   	cout << endl;

   	vector<double> X;
   	X.push_back(-0.25);
   	for (int i = 0; i < RANK - 1 ; i++ ) {
       	X.push_back(0.0);
   	}

   	cout << endl << "X Vector Elements : " << endl;
   	for ( int i = 0 ; i < RANK ; i++ ){
       	cout << X[i] << "   ";
   	}
   	cout << endl;

   	// testing matrix product
   	vector<double> mp;
   	for ( int i = 0 ; i < RANK ; i++ ){
       	mp.push_back(0.0);
   	}
   	compressed::matrixProduct( &mp , &mat_a , &B );
   	cout << endl << "mp Vector Elements : " << endl;
   	for ( int i = 0 ; i < RANK ; i++ ){
       	cout << mp[i] << "   ";
   	}
   	cout << endl << endl;

   	// reset mp to all zeros
   	for ( int i = 0 ; i < RANK ; i++ ){
   		mp[i] = 0.0;
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
           	cout << X[i] << "   ";
       	}
       	cout << " : " << normCurrent << endl;
       	counter++;
   	} 

   	cout<<endl;
   	return 0;
}