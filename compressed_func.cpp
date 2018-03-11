#include "compressed_func.hpp"

void compressed::rowScale( comp_r_mat* A , int i , int j , double a ){
    /* takes row a*row i and adds it to row j in matrix A
     * changes the original matrix; no copy made
     */

    vector<int> row_i_val, row_i_col,row_j_val, row_j_col;
    vector<double>::iterator j_val_start;
    vector<int>::iterator j_col_start;
    j_val_start = A->value.begin() + A->row_p[j];
    j_col_start = A->col_id.begin() + A->row_p[j];
    
    for( int p = A->row_p[i]; p< A->row_p[i+1]; p++){
//        A->value[p] = a* A->value[p];
        row_i_val.push_back(a*A->value[p]);
        row_i_col.push_back(A->col_id[p]);
    }
    
    for( int p = A->row_p[j]; p< A->row_p[j+1]; p++){
        row_j_val.push_back(A->value[p]);
        row_j_col.push_back(A->col_id[p]);
    }
    int size_j = row_j_val.size();
    int row_adj = 0;
    
    for(int i_id = 0; i_id < row_i_col.size(); i_id++){
        bool col_exists = false;
        for(int j_id = 0; j_id < row_j_col.size();j_id++){
            if(row_i_col[i_id] == row_j_col[j_id]){
                row_j_val[j_id] = row_i_val[i_id] + row_j_val[j_id];
                //row_i_col.erase(row_i_col.begin()+i_id);
                col_exists = true;
            }
        }
        if(!col_exists){
            row_adj++;
            for(int j_id = 0; j_id < row_j_col.size();j_id++){
                if(row_i_col[i_id] < row_j_col[j_id]){
                    if(i_id >= row_j_col.size()){
                        row_j_val.insert(row_j_val.end(), row_i_val[i_id]);
                        row_j_col.insert(row_j_col.end(), row_i_col[i_id]);
                        
                    } else{
                        row_j_val.insert(row_j_val.begin()+j_id, row_i_val[i_id]);
                        row_j_col.insert(row_j_col.begin()+j_id, row_i_col[i_id]);
                    }
                    break;
                }else if( j_id == row_j_col.size() - 1){
                    row_j_col.insert(row_j_col.end(), row_i_col[i_id]);
                    row_j_val.insert(row_j_val.end(), row_i_val[i_id]);
                    break;
                }
            }
        }
    }
    A->value.erase(j_val_start, j_val_start + size_j);
    A->col_id.erase(j_col_start, j_col_start + size_j);
    A->value.insert(j_val_start, row_j_val.begin(), row_j_val.end());
    A->col_id.insert(j_col_start, row_j_col.begin(), row_j_col.end());
    for(int p = j+1; p < A->row_p.size(); p++){
        A->row_p[p] = A->row_p[p] + row_adj;
    }
    return;
}

void compressed::rowPermute(comp_r_mat* A, int i, int j){
    /* swaps rows i and j in matrix A
     * changes the original matrix; no copy made
     */
    int small_row, large_row;
    if( i < j ){
        small_row = i;
        large_row = j;
    }else if (j < i){
        small_row = j;
        large_row = i;
    } else {
        return;
    }
    vector<double>::iterator val_start =A->value.begin();
    vector<int>::iterator col_start =A->col_id.begin();
    
    int small_row_count = (A->row_p[small_row + 1] - A->row_p[small_row]);
    int large_row_count = (A->row_p[large_row + 1] - A->row_p[large_row]);
    
    vector<double> row_l_val, row_s_val;
    vector<int> row_l_col, row_s_col;
    
    for( int p = A->row_p[large_row]; p< A->row_p[large_row+1]; p++){
        row_l_val.push_back(A->value[p]) ;
        row_l_col.push_back(A->col_id[p]);
    }
    
    A->value.erase(val_start+A->row_p[large_row], val_start+A->row_p[large_row+1]);
    A->col_id.erase(col_start+A->row_p[large_row], col_start+A->row_p[large_row+1]);
    
    for( int p = A->row_p[small_row]; p< A->row_p[small_row+1]; p++){
        row_s_val.push_back(A->value[p]);
        row_s_col.push_back(A->col_id[p]);
    }
    
    A->value.erase(val_start+A->row_p[small_row], val_start+A->row_p[small_row+1]);
    A->col_id.erase(col_start+A->row_p[small_row], col_start+A->row_p[small_row+1]);
    
    for(int i = small_row+1; i <= large_row; i++){
        A->row_p[i] = A->row_p[i] - small_row_count + large_row_count;
    }
    
    A->value.insert(val_start+A->row_p[small_row], row_l_val.begin(), row_l_val.end());
    A->value.insert(val_start+A->row_p[large_row], row_s_val.begin(), row_s_val.end());
    A->col_id.insert(col_start+A->row_p[small_row], row_l_col.begin(), row_l_col.end());
    A->col_id.insert(col_start+A->row_p[large_row], row_s_col.begin(), row_s_col.end());
    
    return;

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

double compressed::productAx( vector<double>* b, comp_r_mat* A, vector<double>* x ){
    if(A->row_p.size()-1 == x->size()){
        for(int i = 0; i < A->row_p.size()-1; i++){
            double product = 0.0;
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
            (*b)[i]=(product);
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

int compressed::decomposeMatrix( vector<double>* DS , comp_r_mat* LUS , comp_r_mat* AS ){
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
		(*DS).push_back(retrieveElement( AS , i , i ));
		// change the diagonal elements in LUS to 0
		changeElement( LUS , i , i , 0.0 );
	}

	// negate all the elements in LUS
	scalarMultiple( LUS , -1.0 );

	return 0;

}

int compressed::matrixProduct( vector<double>* result , comp_r_mat* A , vector<double>* vec ){
    /* calculates the matrix product of A and vec and stores it in result
     * only works for vector products
     */
    int rank = A->row_p.size() - 1;
    for ( int i = 0 ; i < rank ; i++ ){
        double pdt = 0.0;
		for ( int j = 0 ; j < rank ; j++ ){
			double temp = retrieveElement( A , i , j );
			pdt += (temp)*((*vec)[j]);
		}
        (*result)[i] = pdt;
	}
	return 0;
}

int compressed::jacobiSolver( vector<double>* X , vector<double>* DS , comp_r_mat* LUS , vector<double>* B ){
	int rank = (*X).size();
	vector<double> matPdt;
	for ( int i = 0 ; i < rank ; i ++ )	matPdt.push_back(0.0);

	productAx(&matPdt , LUS , X );
	
	for ( int i = 0 ; i < rank ; i++ ){
		double dInv = 1.0/((*DS)[i]);
		(*X)[i] = ((*B)[i] + matPdt[i])*dInv;
	}

	return 0;

}

int compressed::calculateNorm( double& norm , vector<double>* v , vector<double>* Ax ){
    double squareSum = 0.0;
    int rank = (*v).size();
    for ( int i = 0 ; i < rank ; i++ ){
        double temp = (*v)[i] - (*Ax)[i];
        squareSum +=temp*temp;
    }
    norm = sqrt( squareSum );
    return 0;
}
