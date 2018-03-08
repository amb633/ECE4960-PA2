//
//  main.cpp for programming assignment 2
//
//  Created by Ariana Bruno on 3/7/18.
//  Copyright © 2018 Ariana Bruno. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <fstream>
#include <functional>
#include "compressed_mat_func.hpp"

vector< vector<double> > input_vector = {{1 , 2, 0, 0 , 3}, {4 , 5, 6, 0, 0}, {0 , 7, 8, 0, 9}, {0, 0, 0, 10, 0}, {11 , 0, 0, 0 , 12}};
vector< vector<double> > single_vector = {{5} , {4}, {3}, {2} , {1}};


int main(int argc, const char * argv[]) {
    //Hacker Practice 4.1
    comp_r_mat mat_a = construct_compressed_matrix(&input_vector);
    comp_r_mat mat_b = construct_compressed_matrix(&single_vector);
    print_comp_r_mat( &mat_a );
    comp_r_mat reorder_a = mat_a;
    comp_r_mat reorder_b = mat_b;
    
    // assumes a sqaure matrix
    for( int k = 0; k< mat_a.row_p.size()-1; k++ ){
        reorderMat(&reorder_a, &reorder_a, &reorder_b, k, k);
        print_comp_r_mat(&reorder_a);
        print_comp_r_mat(&reorder_b);
    }

    ifstream rowPtr_file("/Users/arianabruno/Desktop/ECE4960/ProgrammingAssignments/ECE4960-PA2/Mat1/rowPtr.csv");
    int row_value;
    vector<int> row_ptr;
    while(rowPtr_file >> row_value){
        row_ptr.push_back(row_value-1);
    }
    cout << "This is the number of row_ptrs in Mat 1: " << row_ptr.size() << endl;

    ifstream values_file("/Users/arianabruno/Desktop/ECE4960/ProgrammingAssignments/ECE4960-PA2/Mat1/value.csv");
    double value;
    vector<double> values;
    while(values_file >> value){
        values.push_back(value);
    }
    cout << "This is the number of non_zero values in Mat 1: " << values.size() << endl;

    ifstream col_file("/Users/arianabruno/Desktop/ECE4960/ProgrammingAssignments/ECE4960-PA2/Mat1/colInd.csv");
    int col_id;
    vector<int> cols;
    while(col_file >> col_id){
        cols.push_back(col_id-1);
    }
    cout << "This is the number of col_ids for non-zero values in Mat 1: " << cols.size() << endl;
    
    int mat1_rank = row_ptr.size() - 1;
    cout << mat1_rank << endl;
//    comp_r_mat mat1 = construct_compressed_matrix( &row_p, &cols, &values, mat1_rank, mat1_rank);
    comp_r_mat mat1;
    mat1.value = values;
    mat1.col_id = cols;
    mat1.row_p = row_ptr;
    
    
    vector<vector<double>> mat_b_full(mat1_rank);
    for( int i = 0; i< mat_b_full.size(); i++){
        if( i == 1 ) {
            mat_b_full[i] = {1.0};
        } else{
            mat_b_full[i] = {0.0};
        }
    }
    
    comp_r_mat mat_b1 = construct_compressed_matrix(&mat_b_full);
    comp_r_mat reorder_mat1 = mat1;
    comp_r_mat reorder_mat_b1 = mat_b1;
    
    // assumes a sqaure matrix
    for( int k = 0; k< mat1.row_p.size()-1; k++ ){
        reorderMat(&reorder_mat1, &reorder_mat1, &reorder_mat_b1, k, k);
    }
//    print_comp_r_mat(&reorder_mat1);
//    print_comp_r_mat(&reorder_mat_b1);

    return 0;
}
