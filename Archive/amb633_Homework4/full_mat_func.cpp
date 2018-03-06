//
//  full_mat_func.cpp
//  Homework4
//
//  Created by Ariana Bruno on 3/4/18.
//  Copyright © 2018 Ariana Bruno. All rights reserved.
//

#include "full_mat_func.hpp"

double retrieveElementFM( vector<vector<double>>* input, int row_id, int col_id){
    double element = (*input)[row_id][col_id];
    return element;
}

double productAxFM( vector< vector<double>>* input, vector<int>* x, vector<int>* b){
    if((*input)[0].size() == (*x).size()){
        int rows = (*input).size();
        for( int i=0; i<rows; i++ ){
            double product = inner_product((*input)[i].begin(), (*input)[i].end(), (*x).begin(), 0);
            (*b).push_back(product);
        }
        return 0;
    }
    return 1;
}

void mat_scale_FM(vector< vector<double>>* input, int i, int j, double a ){
    vector<double> scale_v;
     
    
    transform ((*input)[i].begin(), (*input)[i].end(), scale_v.begin(), bind1st(multiplies<double>(),a));
    transform (scale_v.begin(), scale_v.end(), (*input)[j].begin(), (*input)[j].begin(), plus<double>());
}

void rowPermute_FM(vector< vector<double>>* input, int i, int j){
    vector<double> hold_v = (*input)[i];
    (*input)[i] = (*input)[j];
    (*input)[j] = hold_v;
}
