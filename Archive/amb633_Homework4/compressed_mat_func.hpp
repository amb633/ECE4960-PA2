//
//  compressed_mat_func.hpp
//  Homework4
//
//  Created by Ariana Bruno on 2/25/18.
//  Copyright © 2018 Ariana Bruno. All rights reserved.
//

#ifndef compressed_mat_func_hpp
#define compressed_mat_func_hpp

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <fstream>

using namespace std;

typedef struct Compressed_Row_Matrix{
    vector<double> value;
    vector<size_t> row_p;
    vector<int> col_id;
}comp_r_mat;


void rowScale(comp_r_mat* A, int i, int j, double a);
void rowPermute(comp_r_mat* A, int i, int j);
double retrieveElement( comp_r_mat* input, int row_id, int col_id);
comp_r_mat construct_compressed_matrix( vector<vector<double>> input );
comp_r_mat construct_compressed_matrix(vector<int>* i, vector<int>* j, vector<double>* val, int rowRank, int colRank);
double productAx( comp_r_mat* A, vector<double>* x, vector<double>* b );
void print_comp_r_mat( comp_r_mat* mat_a );

#endif /* compressed_mat_func_hpp */
