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

void test_compressed::call_tests(){
    cout<<boolalpha;
    
    for (int i = 0; i<5; i++){
        srand(i);
        int rand_i = rand() % test_vector.size();
        int rand_j = rand() % test_vector[0].size();
        cout << "Retrieved correct element " << rand_i << "," << rand_j << " from compressed matrix : " << test_retrieveElement(&test_comp, rand_i, rand_i) << endl;
    }
    
    cout << endl;
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
//            cout << (RS_a*test_vector[RS_i][k])+test_vector[RS_j][k] << " ";
        }
//        cout << endl;
        
        cout << test_rowScale(&scale_test_comp, RS_i, RS_j, RS_a, &correct_row) << endl;
    }
    cout << endl;
    
    for (int i = 0; i<5; i++){
        comp_r_mat row_permute_test_comp = test_comp;
        srand(i);
        int rand_i = rand() % test_vector.size();
        int rand_j = rand() % test_vector[0].size();
        cout << "Row permute accurately switched rows "<< rand_i << " and " << rand_j << ": " << test_rowPermute(&row_permute_test_comp, rand_i, rand_j) << endl;
    }
    
    cout << endl;
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

bool test_compressed::test_construct_compressed_matrix( vector<vector<double>>* input ){
    return true;
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
    return true;
}

bool test_compressed::test_scalarMultiple( comp_r_mat* A , double scale ){
    return true;
}

bool test_compressed::test_copyMatrix( comp_r_mat* C , comp_r_mat* A ){
    return true;
}

bool test_compressed::test_decomposeMatrix( vector<double>* DS , comp_r_mat* LUS , comp_r_mat* AS ){
    return true;
}

bool test_compressed::test_matrixProduct ( vector<double>* result , comp_r_mat* A , vector<double>* vec){
    return true;
}

bool test_compressed::test_jacobiSolver( vector<double>* X , vector<double>* DS , comp_r_mat* LUS , vector<double>* B ){
    return true;
}

bool test_compressed::test_calculateNorm( double& norm , vector<double>* v , vector<double>* Ax ){
    return true;
}

