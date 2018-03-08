#include "sparse.hpp"
#include "full.hpp"

using namespace std;

int main(int argc, char const *argv[])
{
	sparse::sparse_matrix* AS = new sparse::sparse_matrix;
	AS->noofRows = 5;
	AS->noofCols = 5;
	AS->noofVars = 15;

	AS->rowPtr = new int[6];
	AS->colIdx = new int[15];
	AS->value_array = new double[15];

	AS->rowPtr[0] = 0;
	AS->rowPtr[1] = 3;
	AS->rowPtr[2] = 6;
	AS->rowPtr[3] = 9;
	AS->rowPtr[4] = 12;
	AS->rowPtr[5] = 15;

	AS->colIdx[0] = 0;
	AS->colIdx[1] = 1;
	AS->colIdx[2] = 4;
	AS->colIdx[3] = 0;
	AS->colIdx[4] = 1;
	AS->colIdx[5] = 2;
	AS->colIdx[6] = 1;
	AS->colIdx[7] = 2;
	AS->colIdx[8] = 3;
	AS->colIdx[9] = 2;
	AS->colIdx[10] = 3;
	AS->colIdx[11] = 4;
	AS->colIdx[12] = 0;
	AS->colIdx[13] = 3;
	AS->colIdx[14] = 4;

	AS->value_array[0] = -4.0;
	AS->value_array[1] = 1.0;
	AS->value_array[2] = 1.0;
	AS->value_array[3] = 4.0;
	AS->value_array[4] = -4.0;
	AS->value_array[5] = 1.0;
	AS->value_array[6] = 1.0;
	AS->value_array[7] = -4.0;
	AS->value_array[8] = 1.0;
	AS->value_array[9] = 1.0;
	AS->value_array[10] = -4.0;
	AS->value_array[11] = 1.0;
	AS->value_array[12] = 1.0;
	AS->value_array[13] = 1.0;
	AS->value_array[14] = -4.0;

	cout << endl;

	cout << " Original Matrix : " << endl;
	for ( int i = 0 ; i < 5 ; i++ ){
		for ( int j = 0 ; j < 5 ; j++ ){
			double temp;
			sparse::retrieveElement( temp , AS , i , j );
			cout << temp << "   ";
		}
		cout << endl;
	}
	cout << endl;

	sparse::sparse_diagonal_matrix* DS = new sparse::sparse_diagonal_matrix;
	sparse::sparse_matrix* LUS = new sparse::sparse_matrix;

	sparse::decomposeMatrix( DS , LUS , AS );

	cout << " L+U Matrix : " << endl;
	for ( int i = 0 ; i < 5 ; i++ ){
		for ( int j = 0 ; j < 5 ; j++ ){
			double temp;
			sparse::retrieveElement( temp , LUS , i , j );
			cout << temp << "   ";
		}
		cout << endl;
	}
	cout << endl;

	cout << " Diagonal Elements : " << endl;
	for ( int i = 0 ; i < 5 ; i++ ) {
		cout << DS->value_array[i] << "   ";
	}
	cout << endl << endl;

	double* B = new double[5];
	for ( int i = 0 ; i < 5 ; i++ ){
		B[i] = 0.0;
	}
	B[0] = 1.0;

	double* X = new double[5];
	for ( int i = 0 ; i < 5 ; i++ ){
		X[i] = 0.0;
	}
	X[0] = -0.25;

	int counter = 0;
	for ( int j = 0 ; j < 50 ; j++ )
	{
		counter++;
		sparse::jacobiSolver( X , DS , LUS , B );
		cout << counter << " : " ;
		for ( int i = 0 ; i < 5 ; i++ ){
			cout << X[i] << "   ";
		}
		cout << endl;
	}
	return 0;
}