//
//  addtl_funct.cpp
//  program2
//
//  Created by Ariana Bruno on 3/11/18.
//  Copyright © 2018 Ariana Bruno. All rights reserved.
//

#include "addtl_funct.hpp"

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
