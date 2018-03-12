#include "full_func.hpp"


/* returns element at i,j in matrix AF */
double full::retrieveElement( vector< vector<double> >* AF , int rowInd , int colInd ){
	return (*AF)[rowInd][colInd];
}

/* returns element i in vector VF */
double full::retrieveElement( vector<double>* VF , int rowInd ){
	return (*VF)[rowInd];
}

/* print the full matrix AF */
void full::print_full_mat( vector< vector<double>>* AF ){
	// assuming square matrix
	for ( int i = 0 ; i < (*AF).size() ; i++ ){
		for ( int j = 0 ; j < (*AF).size() ; j++ ){
			cout << (*AF)[i][j] << "   ";
		}
		cout << endl;
	}
}

/* print the full vector VF */
void full::print_full_vec( vector<double>* VF ){
	for ( int i = 0 ; i < (*VF).size() ; i++ ){
		cout << (*VF)[i] << "   ";
	}
	cout << endl;
}

/* change an element in matrix AF */
int full::changeElement( vector< vector<double>>* AF , int rowInd , int colInd , double newValue ){
	(*AF)[rowInd][colInd] = newValue;
	return 0;
}

/* create a copy of matrix AF and store it in CF */
int full::copyMatrix ( vector< vector<double>>* CF , vector< vector<double>>* AF ){
	(*CF) = (*AF);
	return 0;
}

/* multiply all (non-zero) elements in AF with scale */
int full::scalarMultiple( vector< vector<double>>* AF , double scale ){
	// assuming square matrix 
	for ( int i = 0 ; i < (*AF).size() ; i++ ){
		for ( int j = 0 ; j < (*AF).size() ; j++ ){
			if ((*AF)[i][j] !=0) (*AF)[i][j] *= scale;
		}
	}
	return 0;
}

/* swap row i,j in matrix AF */
int full::rowPermute( vector< vector<double>>* AF , int i , int j ){
	vector<double> temp = (*AF)[i];
	(*AF)[i] = (*AF)[j];
	(*AF)[j] = temp;
	return 0;
}

/* returns row_j = row_j + a * row_i */
int full::rowScale(vector< vector<double>>* AF, int i, int j, double a ){
	// assuming square matrix
	vector<double> temp = (*AF)[i];
	for ( int ii = 0 ; ii < (*AF).size() ; ii++ ){
		(*AF)[j][ii] += a*temp[ii];
	}
	return 0;

}

/* calculates the matrix product of A and vec and stores it in result */
int full::productAx( vector<double>* result , vector< vector<double>>* AF , vector<double>* VF ){
	// assuming AF is a square matrix, and result and VF are vectors of the same rank as AF
	for ( int i = 0 ; i < (*VF).size() ; i++ ){
		double p = 0.0;
		for ( int j = 0 ; j < (*VF).size() ; j++ ){
			p += ((*AF)[i][j]) * ((*VF)[j]);
		}
		(*result)[i] = p;
	}
	return 0;
}

/* calculate norm of v - AX */
int full::calculateNorm( double& norm , vector<double>* v , vector<double>* Ax ){
	// assuming v and Ax are vectors of the same rank
	double squareSum = 0.0;
	for ( int i = 0 ; i < (*v).size() ; i++ ){
		double temp = (*v)[i] - (*Ax)[i];
		squareSum += temp*temp;
	}
	norm = sqrt( squareSum );
    return 0;
}

/* decompose AF matrix into diagonal elements (stored in DF) and non-diagonal elements (stored in LUF) */
int full::decomposeMatrix( vector<double>* DF , vector< vector<double>>* LUF , vector< vector<double>>* AF ){
	
	copyMatrix( LUF , AF );

	for ( int i = 0 ; i < (*AF).size() ; i++ ){
		// extract the diagonal elements to DS
		(*DF).push_back( retrieveElement( AF , i , i ) );
		// change the diagonal elements in LUF to 0
		changeElement( LUF , i , i , 0.0 );
	}

	// negate all the elements in LUF
	scalarMultiple( LUF , -1.0 );

	return 0;
}

int full::jacobiSolver( vector<double>* X , vector<double>* DF , vector< vector<double>>* LUF , vector<double>* BF ){

	vector<double> matPdt;
	for ( int i = 0 ; i < (*DF).size() ; i++ ) matPdt.push_back(0.0);
	productAx( &matPdt , LUF , X );

	for ( int i = 0 ; i < (*DF).size() ; i++ ){
		double dInv = 1.0/((*DF)[i]);
		(*X)[i] = ((*BF)[i] + matPdt[i])*dInv;
	}
	return 0;
}
