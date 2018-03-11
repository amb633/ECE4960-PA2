#include "compressed_func.hpp"

void compressed::rowScale( comp_r_mat* A , int i , int j , double a ){
    /* takes row a*row i and adds it to row j in matrix A
     * changes the original matrix; no copy made
     */

    /* declares new vectors to hold non-zero elements and corresponding col_id for rows i and j */
    vector<int> row_i_val, row_i_col,row_j_val, row_j_col;
    vector<double>::iterator j_val_start;
    vector<int>::iterator j_col_start;
    j_val_start = A->value.begin() + A->row_p[j];
    j_col_start = A->col_id.begin() + A->row_p[j];
    
    /* multiplies all the i row values by the scalar a */
    for( int p = A->row_p[i]; p< A->row_p[i+1]; p++){
        row_i_val.push_back(a*A->value[p]);
        row_i_col.push_back(A->col_id[p]);
    }
    
    for( int p = A->row_p[j]; p< A->row_p[j+1]; p++){
        row_j_val.push_back(A->value[p]);
        row_j_col.push_back(A->col_id[p]);
    }
    
    int size_j = row_j_val.size();
    int row_adj = 0;
    
    /* checks if the cols in i exist in j */
    /* if col exists in both, the values are added together */
    /* if an i col does NOT exists in j, the value and col_id is inserted into the row j */
    for(int i_id = 0; i_id < row_i_col.size(); i_id++){
        bool col_exists = false;
        for(int j_id = 0; j_id < row_j_col.size();j_id++){
            if(row_i_col[i_id] == row_j_col[j_id]){
                row_j_val[j_id] = row_i_val[i_id] + row_j_val[j_id];
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
    
    /* earses the previous row j values and col_id in matrix A */
    /* inserts in new values and col_id from in row_j_val and row_j_col_id */
    A->value.erase(j_val_start, j_val_start + size_j);
    A->col_id.erase(j_col_start, j_col_start + size_j);
    A->value.insert(j_val_start, row_j_val.begin(), row_j_val.end());
    A->col_id.insert(j_col_start, row_j_col.begin(), row_j_col.end());
    
    /* adjusts row_p for A for the new cols and values added to row j */
    for(int p = j+1; p < A->row_p.size(); p++){
        A->row_p[p] = A->row_p[p] + row_adj;
    }
    return;
}

void compressed::rowPermute(comp_r_mat* A, int i, int j){
    /* swaps rows i and j in matrix A
     * changes the original matrix; no copy made
     */
    
    /* determine which row occurs first in matrix (small_row) */
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
    
    /* create iterator pointer to the beginning of each rows in row_p */
    vector<double>::iterator val_start =A->value.begin();
    vector<int>::iterator col_start =A->col_id.begin();
    int small_row_count = (A->row_p[small_row + 1] - A->row_p[small_row]);
    int large_row_count = (A->row_p[large_row + 1] - A->row_p[large_row]);
    
    /* these vectors will be used to track values and col_id that will be swapped */
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
    
    /* update the row_p to swap the rows, adjust the count in row_p for each row */
    for(int i = small_row+1; i <= large_row; i++){
        A->row_p[i] = A->row_p[i] - small_row_count + large_row_count;
    }
    
    /* insert in the row values and col_id back into the matrix at the swapped starting iterator pointers */
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

    /* finds the non-zero values and push into the A.value vector and the col id into the A.col_id */
    /* counts non-zero values for each row, and pushes the value to row_p for each row (cumulatively) */
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
    /* checks if the dimensions of vector and matrix are the same -- if not return 1*/
    /* checks if product results in inf/ninf -- if yes return 2*/
    /* checks if product results in nan -- if yes return 3*/
    /* checks if product results in sucessful -- if yes return 0*/
    /* product returns to vector b (pre-determined size --same as x) */
    /* product calculated as each row in matrix times the corresponding value in vector and summed */
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

int compressed::jacobiSolver( vector<double>* X , vector<double>* DS , comp_r_mat* LUS , vector<double>* B ){
    int rank = (*X).size();
	vector<double> matPdt;
	for ( int i = 0 ; i < rank ; i ++ )	matPdt.push_back(0.0);
    // calculate product of lower-upper matrix and X
	productAx(&matPdt , LUS , X );
	
    // reverse calculate convergent solution
	for ( int i = 0 ; i < rank ; i++ ){
		double dInv = 1.0/((*DS)[i]);
        // update X values
		(*X)[i] = ((*B)[i] + matPdt[i])*dInv;
	}

	return 0;

}

int compressed::calculateNorm( double& norm , vector<double>* v , vector<double>* Ax ){
    //calculates || v - Ax ||
    double squareSum = 0.0;
    int rank = (*v).size();
    for ( int i = 0 ; i < rank ; i++ ){
        double temp = (*v)[i] - (*Ax)[i];
        squareSum +=temp*temp;
    }
    norm = sqrt( squareSum );
    return 0;
}
