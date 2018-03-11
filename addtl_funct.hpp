//
//  addtl_funct.hpp
//  program2
//
//  Created by Ariana Bruno on 3/11/18.
//  Copyright © 2018 Ariana Bruno. All rights reserved.
//

#ifndef addtl_funct_hpp
#define addtl_funct_hpp

#include <stdio.h>
#include "compressed_func.hpp"

int compressed::matrixProduct( vector<double>* result , comp_r_mat* A , vector<double>* vec );
void compressed::reorderMat( comp_r_mat* input, comp_r_mat* reorder_A, comp_r_mat* reorder_B, int R, int C);
void compressed::columnPermute(comp_r_mat* A, int col1, int col2);
#endif /* addtl_funct_hpp */
