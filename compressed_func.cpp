#include "compressed_func.hpp"

void compressed::rowScale( comp_r_mat* A , int i , int j , double a ){
    /* takes row a*row i and adds it to row j in matrix A
     * changes the original matrix; no copy made 
     */

}

void compressed::rowPermute(comp_r_mat* A, int i, int j){
    /* swaps rows i and j in matrix A 
     * changes the original matrix; no copy made 
     */


}

double compressed::retrieveElement( comp_r_mat* input, int row_id, int col_id){
	/* returns the value stored at specified row_id and col_id */

    double element = 0.0;

    // find cumulative number of elements in previous row
    int row_non_zero_start = input->row_p[row_id];
    // find cumulative number of elements in current row
    int row_non_zero_end = input->row_p[row_id + 1];

    // iterate over current row, until desired col is found
    // if found, return value stored in that index
    for(int p = row_non_zero_start; p< row_non_zero_end; p++){
        if( input->col_id[p] == col_id ){
            element = input->value[p];
            break;
        }
    }
    // if not found, return 0
    return element;
}

compressed::comp_r_mat compressed::construct_compressed_matrix( vector<vector<double>>* input ){
	comp_r_mat mat_A;
	int rows = (*input).size();
	int columns = (*input)[0].size();
	int non_zeros = 0;
	mat_A.row_p.push_back(non_zeros);

	for ( int i = 0 ; i < rows ; i++ ){
		for ( int j = 0 ; j < columns ; j++ ){
			if ((*input)[i][j] !=0){
				mat_A.value.push_back((*input)[i][j]);
				mat_A.col_id.push_back(j);
				non_zeros++;
			}
		}
		mat_A.row_p.push_back(non_zeros);
	}
	return mat_A;
}

compressed::comp_r_mat compressed::construct_compressed_matrix(vector<int>* i, vector<int>* j, vector<double>* val, int rowRank, int colRank){
    comp_r_mat mat_A;
    vector<int> row_p_init(rowRank+1, 0);
    mat_A.row_p = row_p_init;
    
    for( int p = 0; p<(*i).size(); p++){
        int row_p_id = (*i)[p];
        double val_p = (*val)[p];
        int col_p = (*j)[p];
        mat_A.value.insert(mat_A.value.begin() + mat_A.row_p[row_p_id], val_p);
        mat_A.col_id.insert(mat_A.col_id.begin() + mat_A.row_p[row_p_id], (col_p));
        for( int r=row_p_id; r<mat_A.row_p.size(); r++){
            mat_A.row_p[r]++;
        }
    }
    mat_A.row_p.insert(mat_A.row_p.begin(), 0); // offset for the memplus data starting with value 1
    return mat_A;
}

double compressed::productAx( comp_r_mat* A, vector<double>* x, vector<double>* b ){
    if(A->row_p.size()-1 == x->size()){
        for(int i = 0; i < A->row_p.size()-1; i++){
            double product = 0;
            for(int j = 0; j < A->row_p.size()-1; j++){
                double value = (retrieveElement(A, i, j));
                value =  value * (*x)[j];
                if( isinf(value) || isinf(-1*value)){
                    return 2;
                } else if ( isnan(value) ){
                    return 3;
                }
                product = product + value;
            }
            b->push_back(product);
        }
        return 0;
    }
    return 1;

}

void compressed::print_comp_r_mat( comp_r_mat* mat_a ){
    cout << "This is the values in mat: " << endl;
    for (double n = 0; n<mat_a->value.size() ; n++){
        cout << mat_a->value[n] << ", ";
    }
    cout << endl;
    cout << "This is the row_p in mat: " << endl;
    for (double n = 0; n<mat_a->row_p.size() ; n++){
        cout << mat_a->row_p[n] << ", ";
    }
    cout << endl;
    cout << "This is the col_id in mat: " << endl;
    for (double n = 0; n<mat_a->col_id.size() ; n++){
        cout << mat_a->col_id[n] << ", ";
    }
    cout << endl;

}

bool compressed::check_sum( comp_r_mat* mat, vector<double>* vec ){
    bool check_sum = false;
    
    double mat_sum = 0;
    for(int p=0; p<mat->value.size(); p++){
        mat_sum = mat_sum + mat->value[p];
    }
    double vec_sum = 0;
    for(int k=0; k<(*vec).size(); k++){
        vec_sum = vec_sum + (*vec)[k];
    }
    double sum_dif = abs(mat_sum-vec_sum);
    if(sum_dif <= (pow(10.0,-7))){
        cout << sum_dif << endl;
        check_sum=true;
    }
    return check_sum;

}

void compressed::reorderMat( comp_r_mat* input, comp_r_mat* reorder_A, comp_r_mat* reorder_B, int R, int C){
    int row_id = R;
    double max = 0.0;
    int max_R = 0;
    int max_C = 0;
    
    int v_counter = 0;
    vector<int>* col_id = &(input->col_id);
    vector<double>* value = &(input->value);
    int row_start = input->row_p[R];
    
    for( int i = row_start; i<col_id->size(); i++){
        int row_non_zero_end = input->row_p[row_id + 1];
        if(i >= row_non_zero_end){
            row_id++;
        }
        if((*col_id)[i] >= C){
            if((*value)[i] > max){
                max = (*value)[i];
                v_counter = i;
                max_R = row_id;
                max_C = (*col_id)[i];
            }
        }
    }
    if (R != max_R){
        cout << "This is the next max value: " << max << " and row swap " << R << " with " << max_R << endl;
        rowPermute(reorder_A, max_R, R);
        rowPermute(reorder_B, max_R, R);
    }
    if(C != max_C){
        cout << "This is the next max value: " << max << " and col swap " << C << " with " << max_C << endl;
        columnPermute(reorder_A, max_C, C);
    }

}

void compressed::columnPermute(comp_r_mat* A, int col1, int col2){
    for( int i= 1; i< A->row_p.size(); i++){
        int row_values = A->row_p[i];
        int swap1_id = -1, swap2_id = -1;
        for(int j = A->row_p[i-1]; j < row_values; j++ ){
            if(A->col_id[j] == col1){
                swap1_id = j;
            }else if(A->col_id[j] == col2){
                swap2_id = j;
            }
        }
        if( swap1_id >= 0 && swap2_id >= 0 ){
            double hold_val = A->value[swap1_id];
            
            A->value[swap1_id] = A->value[swap2_id];
            
            A->value[swap2_id] = hold_val;
        } else if( swap1_id >= 0 || swap2_id >= 0 ){
            int non_zero_id = swap1_id;
            int zero_id = swap2_id;
            int col = col1;
            int zero_col = col2;
            if(swap2_id >=0 ){
                non_zero_id = swap2_id;
                zero_id = swap1_id;
                col = col2;
                zero_col = col1;
            }
            double value = A->value[non_zero_id];
            for(int j = A->row_p[i-1]; j < row_values; j++ ){
                if(A->col_id[j] >= zero_col || j == (row_values-1)){
                    (A->col_id).erase(A->col_id.begin() + non_zero_id);
                    (A->value).erase(A->value.begin() + non_zero_id);
                    (A->col_id).insert(A->col_id.begin()+j, zero_col);
                    (A->value).insert(A->value.begin()+j, value);
                    break;
                }
                
            }
        }
    }
}

int compressed::changeElement( comp_r_mat* A , int rowInd , int colInd , double newValue ){
	/* changes the value of an existing non-zero element
	 * returns 0 if the element was successfully changed 
	 * returns -1 if the element is currently a zero, and couldn't be changed 
	 */
	int start = A->row_p[rowInd];
	int end = A->row_p[rowInd+1];
	for ( int i = start ; i < end ; i++ ){
		if ( A->col_id[i] == colInd ){
			A->value[i] = newValue;
			return 0;
		}
	}
	return -1;
}

int compressed::scalarMultiple( comp_r_mat* A , double scale ){
	/*  scalar operations (e.g. want to multiply all elements in matrix by 2);
	 *  so the row and column order are not important, just muliply everything
	 */
	for ( int i = 0 ; i < A->value.size() ; i++ ){
		if( A->value[i] != 0 )	A->value[i] = scale*A->value[i];
	}
	return 0;
}

int compressed::copyMatrix( comp_r_mat* C , comp_r_mat* A ){
	/* creates a copy of matrix A in matrix C
	 * create a pointer to both in the main, and pass the pointers as function handles
	 * need to do a deep copy if using dynamic memory allocation
	 */
	/*C->noofRows = A->row_p.size() - 1;
	C->noofCols = A->row_p.size() - 1;
	C->noofVars = A->value.size();*/
	C->value = A->value;
	C->row_p = A->row_p;
	C->col_id = A->col_id;
	return 0;

}

int compressed::decomposeMatrix( comp_diag* DS , comp_r_mat* LUS , comp_r_mat* AS ){
	/* decomposes matrix AS into its diagonal elements (stored in DS)
	 * and its lower upper form (stored in LUS)
	 * create all three structures in main, and pass in as function handles
	 */

	// create a copy of AS in LUS
	copyMatrix( LUS , AS );

	// DS is a single vector of doubles containing the diagonal elements
    int rank = AS->row_p.size() - 1;

	for ( int i = 0 ; i < rank ; i++ ){
		// extract the diagonal elements to DS
		DS->value.push_back(retrieveElement( AS , i , i ));
		// change the diagonal elements in LUS to 0
		changeElement( LUS , i , i , 0.0 );
	}

	// negate all the elements in LUS
	scalarMultiple( LUS , -1.0 );

	return 0;

}

int compressed::matrixProduct( comp_diag* result , comp_r_mat* A , comp_diag* vec ){
    /* calculates the matrix product of A and vec and stores it in result
     * only works for vector products
     */
    int rank = A->row_p.size() - 1;
    for ( int i = 0 ; i < rank ; i++ ){
        double pdt = 0.0;
		for ( int j = 0 ; j < rank ; j++ ){
			double temp = retrieveElement( A , i , j );
			pdt += temp*(vec->value[j]);
		}
        result->value[i] = pdt;
	}
	return 0;
}

int compressed::jacobiSolver( comp_diag* X , comp_diag* DS , comp_r_mat* LUS , comp_diag* B ){
	int rank = X->value.size();
	comp_diag matPdt;
	for ( int i = 0 ; i < rank ; i ++ )	matPdt.value.push_back(0.0);

	matrixProduct(&matPdt , LUS , X );
	
	for ( int i = 0 ; i < rank ; i++ ){
		double dInv = 1.0/DS->value[i];
		X->value[i] = (B->value[i] + matPdt.value[i])*dInv;
	}

	return 0;

}

int compressed::calculateNorm( double& norm , comp_diag* v , comp_diag* Ax ){
    double squareSum = 0.0;
    for ( int i = 0 ; i < v->value.size() ; i++ ){
        double temp = v->value[i] - Ax->value[i];
        squareSum +=temp*temp;
    }
    norm = sqrt( squareSum );
    return 0;
}
