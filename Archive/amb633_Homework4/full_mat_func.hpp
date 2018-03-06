//
//  full_mat_func.hpp
//  Homework4
//
//  Created by Ariana Bruno on 3/4/18.
//  Copyright © 2018 Ariana Bruno. All rights reserved.
//

#ifndef full_mat_func_hpp
#define full_mat_func_hpp

#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <stdio.h>

using namespace std;

double retrieveElementFM( vector<vector<double>>* input, int row_id, int col_id);
double productAxFM( vector< vector<double>>* input, vector<int>* x, vector<int>* b);
void mat_scale_FM(vector< vector<double>>* input, int i, int j, double a );
void rowPermute_FM(vector< vector<double>>* input, int i, int j);

#endif /* full_mat_func_hpp */
