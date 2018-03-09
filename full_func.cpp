#include "full_func.hpp"

double full::retrieveElement( vector< vector<double> >* A , int row_id , int col_id ){
	double element = (*A)[row_id][col_id];
	return element;
}