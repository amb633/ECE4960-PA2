//
//  test_comp_func.cpp
//  program2
//
//  Created by Ariana Bruno on 3/9/18.
//  Copyright © 2018 Ariana Bruno. All rights reserved.
//

#include "test_comp_func.hpp"
#include <stdlib.h>

vector<vector<double>> test_vector = {{-4, 1, 0, 0, 1}, {4, -4, 1, 0, 0}, { 0, 1, -4, 1, 0}, {0, 0, 1, -4, 1}, {1, 0, 0, 1, -4}};
comp_r_mat test_comp = construct_compressed_matrix(&test_vector);

//using namespace compressed;

void test_compressed::call_comp_tests(){
    cout<<boolalpha;
    
    /* Testing the retrieveElement function on random indices of the matrix (multiple iterations) */
    for (int i = 0; i<5; i++){
        srand(i);
        int rand_i = rand() % test_vector.size();
        int rand_j = rand() % test_vector[0].size();
        cout << "Retrieved correct element " << rand_i << "," << rand_j << " from compressed matrix : " << test_retrieveElement(&test_comp, rand_i, rand_i) << endl;
    }
    cout << endl;
    
    /* Testing the rowScale function on random rows of the matrix (multiple iterations) */
    for (int i = 0; i<5; i++){
        comp_r_mat scale_test_comp = test_comp;
        srand(i);
        int RS_i = rand() % test_vector.size();
        int RS_j = rand() % test_vector[0].size();
        double RS_a = (rand() % 4)*1.0;

        cout << "Testing compressed matrix rowScale: row " << RS_i << " multiplied by " << RS_a << " and added to row " << RS_j << " returns correct result: ";
        
        vector<double> correct_row = {};
        for( int k =0; k<test_vector[0].size(); k++){
            correct_row.push_back((RS_a*test_vector[RS_i][k])+test_vector[RS_j][k]);
        }
        
        cout << test_rowScale(&scale_test_comp, RS_i, RS_j, RS_a, &correct_row) << endl;
    }
    cout << endl;
    
    /* Testing the rowPermute function on random rows of the matrix (multiple iterations) */
    for (int i = 0; i<5; i++){
        comp_r_mat row_permute_test_comp = test_comp;
        srand(i);
        int rand_i = rand() % test_vector.size();
        int rand_j = rand() % test_vector[0].size();
        cout << "Row permute accurately switched rows "<< rand_i << " and " << rand_j << ": " << test_rowPermute(&row_permute_test_comp, rand_i, rand_j) << endl;
    }
    cout << endl;
    
    /* Testing the productAx function using given vectors and expected result vectors (multiple iterations) */
    vector<double> x_0 = {0, 0, 0, 0, 0};
    vector<double> expected_b_0 = {0, 0, 0, 0, 0};
    cout << "Multiplying matrix by vector: { ";
    for(int i = 0; i<x_0.size(); i++){
        cout << x_0[i] << " ";
    }
    cout << "} results in a vector: { ";
    for(int i = 0; i<x_0.size(); i++){
        cout << expected_b_0[i] << " ";
    }
    cout << "} = " << test_productAx(&test_comp, &x_0, &expected_b_0) << endl;
    
    
    vector<double> x_1 = {1, 0, 0, 0, 0};
    vector<double> expected_b_1 = {-4, 4, 0, 0, 1};
    cout << "Multiplying matrix by vector: { ";
    for(int i = 0; i<x_1.size(); i++){
        cout << x_1[i] << " ";
    }
    cout << "} results in a vector: { ";
    for(int i = 0; i<expected_b_1.size(); i++){
        cout << expected_b_1[i] << " ";
    }
    cout << "} = " << test_productAx(&test_comp, &x_1, &expected_b_1) << endl;
    
    vector<double> x_2 = {1, 0, 0, 0, 1};
    vector<double> expected_b_2 = {-3, 4, 0, 1, -3};
    cout << "Multiplying matrix by vector: { ";
    for(int i = 0; i<x_2.size(); i++){
        cout << x_2[i] << " ";
    }
    cout << "} results in a vector: { ";
    for(int i = 0; i<expected_b_2.size(); i++){
        cout << expected_b_2[i] << " ";
    }
    cout << "} = " << test_productAx(&test_comp, &x_0, &expected_b_0) << endl;
    cout << endl;
    
    /* Testing the changeElement function on random indicies of the matrix (multiple iterations) */
    for (int i = 0; i<5; i++){
        comp_r_mat change_element_test_comp = test_comp;
        srand(i);
        int rand_i = rand() % test_vector.size();
        int rand_j = rand() % test_vector[0].size();
        double rand_val = (rand() % 10)*1.0;
        cout << "Change element accurately switched updates element at ("<< rand_i << "," << rand_j << "): " << test_changeElement(&change_element_test_comp, rand_i, rand_j, rand_val)<< endl;
    }
    cout << endl;
    
    /* Testing the scalarMultiple function on random scalar values to scale the matrix (multiple iterations) */
    for (int i = 0; i<5; i++){
        comp_r_mat scale_test_comp = test_comp;
        srand(i);
        double rand_val = (rand() % 10)*1.0;
        cout << "Accurately scale matrix by " << rand_val << ": " << test_scalarMultiple(&scale_test_comp, rand_val)<< endl;
    }
    cout << endl;
    
    /* Testing the copyMatrix function on matrix (single test) */
    comp_r_mat copy_test_comp;
    cout << "Accurately copy matrix: " << test_copyMatrix(&copy_test_comp, &test_comp);
    cout << endl;
    
    /* Testing the decomposeMatrix function on matrix (single test for correct and two purposely failed tests) */
    cout << endl;
    vector<vector<double>> expected_LUS_vec = {{0, -1, 0, 0, -1}, {-4, 0, -1, 0, 0}, { 0, -1, 0, -1, 0}, {0, 0, -1, 0, -1}, {-1, 0, 0, -1, 0}};
    comp_r_mat expected_LUS_comp = construct_compressed_matrix(&expected_LUS_vec);
    vector<double> expected_DS = {-4, -4, -4, -4, -4};
    cout << "Accurately decompose matrix: " << test_decomposeMatrix(&expected_DS, &expected_LUS_comp, &test_comp);
    cout << endl;
    vector<double> incorrect_DS = {-1, -4, -4, -4, -4};
    cout << "Test incorrect decompose DS matrix: " << test_decomposeMatrix(&incorrect_DS, &expected_LUS_comp, &test_comp);
    cout << endl;
    comp_r_mat incorrect_LUS = expected_LUS_comp;
    incorrect_LUS.value[0] = 2*incorrect_LUS.value[0];
    cout << "Test incorrect decompose LUS matrix: " << test_decomposeMatrix(&expected_DS, &incorrect_LUS, &test_comp);
    cout << endl;
    
    /* Testing the copyMatrix function on matrix (single test) */
    cout << endl;
    vector<double> v1_1 = {4.0,2.0,2.0,2.0,2.0};
    vector<double> v1_2 = {4.0,0.0,0.0,0.0,0.0};
    double expected_norm1 = 4.0;
    cout << "Accurately calculates the residual norm between two vectors: " << test_calculateNorm(expected_norm1, &v1_1, &v1_2);
    cout << endl;
    vector<double> v2_1 = {1.0,-1.0,0.0,-1.0,1.0};
    vector<double> v2_2 = {0.0,0.0,0.0,0.0,0.0};
    double expected_norm2 = 2.0;
    cout << "Accurately calculates the residual norm between two vectors: " << test_calculateNorm(expected_norm2, &v2_1, &v2_2);
    cout << endl;
    
}
    

bool test_compressed::test_rowScale(comp_r_mat* input, int i, int j, double a, vector<double>* true_result){
    rowScale(input, i, j, a);
    bool correct = true;
    for( int k = 0; k<true_result->size(); k++){
        double comp_value = retrieveElement(input, j, k);
        if(comp_value != (*true_result)[k]){
            correct = false;
        }
    }
    return correct;
}
    
bool test_compressed::test_rowPermute(comp_r_mat* input, int i, int j){
    rowPermute(input, i, j);
    bool check = true;
    double element;
    for(int k = 0; k<input->row_p.size()-1 ; k++){
        element = retrieveElement(input, i, k);
        if( element != test_vector[j][k]){
            check = false;
        }
    }
    for(int k = 0; k<input->row_p.size()-1 ; k++){
        element = retrieveElement(input, j, k);
        if( element != test_vector[i][k]){
            check = false;
        }
    }
    
    return check;
}

bool test_compressed::test_retrieveElement(comp_r_mat* input, int i, int j){
    /* check test_vector[0][0] returns the correct element in the compressed matrix format */
    double correct_value = test_vector[i][j];
    double comp_value = retrieveElement(input, i, j);
    if (correct_value == comp_value){
        return true;
    }
    return false;
}

bool test_compressed::test_productAx( comp_r_mat* A, vector<double>* x, vector<double>* expected_result ){
    vector<double> b(expected_result->size());
    productAx(&b, A, x);
    bool check = true;
    for( int i = 0; i<b.size(); i++){
        if(b[i] != (*expected_result)[i]){
            check = false;
        }
    }
    return check;
}


bool test_compressed::test_changeElement( comp_r_mat* A , int rowInd , int colInd , double newValue ){
    int error = changeElement(A, rowInd, colInd, newValue);
    if (error == 0){
        double value = retrieveElement(A, rowInd, colInd);
        if(value == newValue){
            return true;
        }
        return false;
    }else{
        cout << "cannot change zero valued elements -- ";
        return true;
    }
}

bool test_compressed::test_scalarMultiple( comp_r_mat* A , double scale ){
    comp_r_mat A2 = *A;
    scalarMultiple(&A2, scale);
    bool check = true;
    for( int i = 0; i<(A2.value.size()); i++){
        if( A2.value[i] != (scale*(A->value[i])) ){
            check = false;
        }
    }
    return check;
}

bool test_compressed::test_copyMatrix( comp_r_mat* C , comp_r_mat* A ){
    copyMatrix(C, A);
    bool check = true;
    for( int i=0; i<C->row_p.size()-1; i++){
        for( int j=0; j< C->row_p.size()-1; j++){
            double c_element = retrieveElement(C, i, j);
            double a_element = retrieveElement(A, i, j);
            if( c_element != a_element ){
                check = false;
            }
        }
    }
    return check;
}

bool test_compressed::test_decomposeMatrix( vector<double>* expected_DS , comp_r_mat* expected_LUS , comp_r_mat* AS ){
    vector<double> DS;
    comp_r_mat LUS;
    decomposeMatrix(&DS, &LUS, AS);
    bool LUS_check = true;
    for( int i=0; i<LUS.row_p.size()-1; i++){
        for( int j=0; j< LUS.row_p.size()-1; j++){
            double result = retrieveElement(&LUS, i, j);
            double expected = retrieveElement(expected_LUS, i, j);
            if( result != expected ){
                LUS_check = false;
            }
        }
    }
    bool DS_check = true;
    for( int k=0; k<DS.size(); k++){
        double result = DS[k];
        double expected = (*expected_DS)[k];
        if( result != expected ){
            DS_check = false;
        }
    }
    
    if(LUS_check && DS_check){
        return true;
    }else{
        if(LUS_check){
            cout << "-- diagonal result are incorrect ";
            return false;
        }else if(DS_check){
            cout << "-- lower/upper results are incorrect ";
            return false;
        }else{
            cout<< "all results are incorrect";
            return false;
        }
            
    }
    return false;
}

bool test_compressed::test_calculateNorm( double& expected_norm, vector<double>* v , vector<double>* Ax ){
    double norm;
    calculateNorm(norm, v, Ax);
    if( norm != expected_norm ){
        return false;
    }
    return true;
}

