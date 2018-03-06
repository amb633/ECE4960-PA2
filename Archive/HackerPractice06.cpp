#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

int jacobiSolver( double* X , double** D , double ** L , double* B );
int matrixProduct( double* b , double** A , double* v );
int calculateNorm( double& norm , double* v , double* Ax );

int main(int argc, char const *argv[])
{
	double** A;
	A = new double*[5];
	for ( int i = 0 ; i < 5 ; i++ ){
		A[i] = new double[5];
	}
	for ( int i = 0 ; i < 5 ; i++ ){
		for ( int j = 0 ; j < 5 ; j++ ){
			A[i][j] = 0.0;
			if ( i == j ) {
				A[i][j] = -4.0;
			}
		}
	}

	A[0][1] = 1.0;
	A[1][2] = 1.0;
	A[2][3] = 1.0;
	A[3][4] = 1.0;
	A[1][0] = 4.0;
	A[2][1] = 1.0;
	A[3][2] = 1.0;
	A[4][3] = 1.0;
	A[4][0] = 1.0;
	A[0][4] = 1.0;

	double** D;
	D = new double*[5];
	for ( int i = 0 ; i < 5 ; i++ ){
		D[i] = new double[5];
	}
	for ( int i = 0 ; i < 5 ; i++ ){
		for ( int j = 0 ; j < 5 ; j++ ){
			D[i][j] = 0.0;
			if ( i == j ) {
				D[i][j] = -4.0;
			}
		}
	}


	double** L;
	L = new double*[5];
	for ( int i = 0 ; i < 5 ; i++ ){
		L[i] = new double[5];
	}
	for ( int i = 0 ; i < 5 ; i++ ){
		for ( int j = 0 ; j < 5 ; j++ ){
			L[i][j] = 0.0;
		}
	}

	L[1][0] = -4.0;
	L[2][1] = -1.0;
	L[3][2] = -1.0;
	L[4][3] = -1.0;
	L[4][0] = -1.0;
	L[0][1] = -1.0;
	L[1][2] = -1.0;
	L[2][3] = -1.0;
	L[3][4] = -1.0;
	L[0][4] = -1.0;

	double* B;
	B = new double[5];
	for( int i = 0 ; i < 5 ; i ++ ){
		B[i] = 0.0;
	}
	B[0] = 1.0;

	double* X;
	X = new double[5];
	for ( int i = 0 ; i < 5 ; i++ ){
		X[i] = 0.0;
	}
	X[0] = -0.25;

	double normCurrent , normPrev;
	normPrev = 0.1;
	double res[5] = { 0.0 , 0.0 , 0.0 , 0.0 , 0.0 };
	double zero[5] = { 0.0 , 0.0 , 0.0 , 0.0 , 0.0 };
	jacobiSolver( X , D , L , B );
	matrixProduct( res , A , X );
	calculateNorm( normCurrent , B , res );
	cout << endl << normCurrent << endl << endl;

	int counter = 0;

	//for( int j = 0 ; j < 50 ; j++ )
	while( abs(normCurrent-normPrev) > 10e-7 ) 
	{
		normPrev = normCurrent;
		counter++;
		jacobiSolver( X , D , L , B );

		double resCurrent[5] = { 0.0 , 0.0 , 0.0 , 0.0 , 0.0 };
		matrixProduct( resCurrent , A , X );
		calculateNorm( normCurrent , B , resCurrent );

		cout << counter << " : " ;
		for ( int i = 0 ; i < 5 ; i++ ){
			cout << X[i] << "   " ;
		}
		cout << " : " << normCurrent << endl;
	}

	cout << endl;
	return 0;
}

int jacobiSolver( double* X , double** D , double ** L , double* B ){
	double dInv[5] = {0.0 , 0.0 , 0.0 , 0.0 , 0.0 };
	double matPdt[5] = {0.0 , 0.0 , 0.0 , 0.0 , 0.0 };
	for ( int i = 0 ; i < 5 ; i++ ){
		dInv[i] = 1.0/D[i][i];
		//cout << dInv[i] << endl;
	}
	matrixProduct( matPdt , L , X );

	for ( int i = 0 ; i < 5 ; i++ ){
		X[i] = (B[i] + matPdt[i])*dInv[i];
		//cout << X[i] << endl;
	}

	return 0;

}

int matrixProduct( double* b , double** A , double* v){
	for ( int i = 0 ; i < 5 ; i++ ){
		for ( int j = 0 ; j < 5 ; j++ ){
			b[i] += A[i][j]*v[j];
		}
	}
}

int calculateNorm( double& norm , double* v , double* Ax ) {
	double squareSum = 0.0;
	for ( int i = 0 ; i < 5 ; i++ ){
		double temp = v[i] - Ax[i];
		squareSum += temp*temp;
	}
	norm = sqrt( squareSum );	
}