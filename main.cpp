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
    
    compressed::comp_r_mat mat1;
    mat1.value = values;
    mat1.col_id = cols;
    mat1.row_p = row_ptr;
    
//    vector<vector<double>> mat_b1_full(mat1_rank);
//    for( int i = 0; i< mat_b1_full.size(); i++){
//        if( i == 0 ) {
//            mat_b1_full[i] = {1.0};
//        } else{
//            mat_b1_full[i] = {0.0};
//        }
//    }
//
//    compressed::comp_r_mat mat_b1_comp = compressed::construct_compressed_matrix(&mat_b1_full);
    
    //    TODO: currently commented out reordering
//    compressed::comp_r_mat reorder_mat1 = mat1;
//    compressed::comp_r_mat reorder_mat_b = mat_b_comp;
//
////    assumes a sqaure matrix
//    for( int k = 0; k< mat1.row_p.size()-1; k++ ){
//        reorderMat(&reorder_mat1, &reorder_mat1, &reorder_mat_b, k, k);
//    }
    
    vector<double> mat_b1;
    for( int i = 0; i< mat1_rank; i++){
        if(i == 0){
            mat_b1.push_back(1.0);
        } else {
            mat_b1.push_back(0.0);
        }
    }
    
    vector<double> mat_b2;
    for( int i = 0; i< mat1_rank; i++){
        if(i == 4){
            mat_b2.push_back(1.0);
        } else {
            mat_b2.push_back(0.0);
        }
    }
    
    vector<double> mat_b3;
    for( int i = 0; i< mat1_rank; i++){
        mat_b3.push_back(1.0);
    }

    
    compressed::comp_r_mat LUC;
    vector<double> DC;
    compressed::decomposeMatrix( &DC , &LUC , &mat1 );
    
    vector<double> mat_x1;
    for( int i = 0; i< mat1_rank; i++){
        mat_x1.push_back(0.0);
    }
    mat_x1[0] = (1.0/DC[0]);
    
    vector<double> matProd1;
    for ( int i = 0 ; i < mat1_rank ; i++ ){
    	matProd1.push_back(0.0);
    }
    
    double normPrev1 = 2;
    double normCurrent1 = 1;
    int counter1 = 0;

    // calculating mat_x for mat_b1
    cout << "calculating mat_x for mat_b1" << endl;
    while( abs(normCurrent1 - normPrev1) > 1e-10 ){
    	normPrev1 = normCurrent1;
        compressed::jacobiSolver( &mat_x1 , &DC , &LUC , &mat_b1 );
        compressed::matrixProduct( &matProd1 , &mat1 , &mat_x1 );
        compressed::calculateNorm( normCurrent1 , &mat_b1 , &matProd1 );
        cout << counter1 << " : " ;
        for ( int i = 0 ; i < 5 ; i++ ){
            cout << mat_x1[i] << "   ";
        }
        cout << " : " << normCurrent1 << endl;
        counter1++;
    }
    
    vector<double> matProd2;
    for ( int i = 0 ; i < mat1_rank ; i++ ){
        matProd2.push_back(0.0);
    }
    
    double normPrev2 = 2;
    double normCurrent2 = 1;
    int counter2 = 0;
    
    vector<double> mat_x2;
    for( int i = 0; i< mat1_rank; i++){
        mat_x2.push_back(0.0);
    }
    mat_x2[4] = (1.0/DC[4]);
    
    // calculating mat_x for mat_b2
    cout << "calculating mat_x for mat_b2" << endl;
    while( abs(normCurrent2 - normPrev2) > 1e-10 ){
        normPrev2 = normCurrent2;
        compressed::jacobiSolver( &mat_x2 , &DC , &LUC , &mat_b2 );
        compressed::matrixProduct( &matProd2 , &mat1 , &mat_x2 );
        compressed::calculateNorm( normCurrent2 , &mat_b2 , &matProd2 );
        cout << counter2 << " : " ;
        for ( int i = 0 ; i < 5 ; i++ ){
            cout << mat_x2[i] << "   ";
        }
        cout << " : " << normCurrent2 << endl;
        counter2++;
    }
    
    vector<double> matProd3;
    for ( int i = 0 ; i < mat1_rank ; i++ ){
        matProd3.push_back(0.0);
    }
    
    double normPrev3 = 2;
    double normCurrent3 = 1;
    int counter3 = 0;

    vector<double> zeros;
    for ( int i = 0 ; i < mat1_rank ; i++ ){
    	zeros.push_back(0.0);
    }
    
    vector<double> mat_x3;
    for( int i = 0; i< mat1_rank; i++){
        mat_x3.push_back((1.0/DC[i]));
    }

    double norm_b;
    compressed::calculateNorm( norm_b , &mat_b3 , &zeros );
    
    cout << "calculating mat_x for mat_b3" << endl;
    while( abs(normCurrent3 - normPrev3) > 1e-10 ){
        normPrev3 = normCurrent3;
        compressed::jacobiSolver( &mat_x3 , &DC , &LUC , &mat_b3 );
        compressed::matrixProduct( &matProd3 , &mat1 , &mat_x3 );
        compressed::calculateNorm( normCurrent3 , &mat_b3 , &matProd3 );
        cout << counter3 << " : " ;
        for ( int i = 0 ; i < 5 ; i++ ){
            cout << mat_x3[i] << "   ";
        }
        cout << " : " << normCurrent3/norm_b << endl;
        counter3++;
    }
    
}
