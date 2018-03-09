//
//  main.cpp
//  
//
//  Created by Ariana Bruno on 3/8/18.
//

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>

#include "compressed_func.hpp"
#include "full_func.hpp"

using namespace std;

int main(int argc, char const *argv[])
{
    ifstream rowPtr_file("./Mat1/rowPtr.csv");
    int row_value;
    vector<int> row_ptr;
    while(rowPtr_file >> row_value){
        row_ptr.push_back(row_value-1);
    }
    cout << "This is the number of row_ptrs in Mat 1: " << row_ptr.size() << endl;
    
    ifstream values_file("./Mat1/value.csv");
    double value;
    vector<double> values;
    while(values_file >> value){
        values.push_back(value);
    }
    cout << "This is the number of non_zero values in Mat 1: " << values.size() << endl;
    
    ifstream col_file("./Mat1/colInd.csv");
    int col_id;
    vector<int> cols;
    while(col_file >> col_id){
        cols.push_back(col_id-1);
    }
    cout << "This is the number of col_ids for non-zero values in Mat 1: " << cols.size() << endl;
    
    int mat1_rank = row_ptr.size() - 1;
    cout << mat1_rank << endl;
    //    comp_r_mat mat1 = construct_compressed_matrix( &row_p, &cols, &values, mat1_rank, mat1_rank);
    compressed::comp_r_mat mat1;
    mat1.value = values;
    mat1.col_id = cols;
    mat1.row_p = row_ptr;
    mat1.noofRows = mat1.row_p.size()-1;
    mat1.noofCols = mat1.row_p.size()-1;
    mat1.noofVars = mat1.value.size();
    
    vector<vector<double>> mat_b_full(mat1_rank);
    for( int i = 0; i< mat_b_full.size(); i++){
        if( i == 0 ) {
            mat_b_full[i] = {1.0};
        } else{
            mat_b_full[i] = {0.0};
        }
    }
    
    compressed::comp_r_mat mat_b_comp = compressed::construct_compressed_matrix(&mat_b_full);
    
    compressed::comp_r_mat reorder_mat1 = mat1;
    compressed::comp_r_mat reorder_mat_b = mat_b_comp;
    
    // assumes a sqaure matrix
    for( int k = 0; k< mat1.row_p.size()-1; k++ ){
        reorderMat(&reorder_mat1, &reorder_mat1, &reorder_mat_b, k, k);
    }
    
    compressed::comp_diag mat_b;
    mat_b.rank = mat1_rank;
    for( int i = 0; i< mat_b.rank; i++){
        mat_b.value.push_back(compressed::retrieveElement(&reorder_mat_b, i, 0));
    }
    
    compressed::comp_diag mat_x;
    mat_x.rank = mat1_rank;
    for( int i = 0; i< mat_x.rank; i++){
        if( i == 0 ) {
            mat_x.value.push_back(-0.25);
        } else{
            mat_x.value.push_back(0.0);
        }
    }

//    compressed::comp_r_mat mat_b2;
    
    compressed::comp_r_mat LUC;
    compressed::comp_diag DC;
    compressed::decomposeMatrix( &DC , &LUC , &mat1 );
    
    int counter = 0;
    for ( int j = 0 ; j < 50 ; j++ )
    {
        counter++;
        compressed::jacobiSolver( &mat_x , &DC , &LUC , &mat_b );
        cout << counter << " : " ;
        for ( int i = 0 ; i < 5 ; i++ ){
            cout << mat_x.value[i] << "   ";
        }
        cout << endl;
    }
    
    
}
