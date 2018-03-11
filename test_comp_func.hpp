//
//  test_comp_func.hpp
//  program2
//
//  Created by Ariana Bruno on 3/9/18.
//  Copyright © 2018 Ariana Bruno. All rights reserved.
//

#ifndef test_comp_func_hpp
#define test_comp_func_hpp

#include <stdio.h>
#include "compressed_func.hpp"

using namespace compressed;

namespace test_compressed{
    
    void call_tests();
    
    bool test_rowScale(comp_r_mat* input, int i, int j, double a, vector<double>* true_result);
    
    bool test_rowPermute(comp_r_mat* input, int i, int j);
    
    bool test_retrieveElement(comp_r_mat* input, int i, int j);
    
    bool test_construct_compressed_matrix( vector<vector<double>>* input );
    
    
    bool test_productAx( comp_r_mat* A, vector<double>* x, vector<double>* expected_result );
    
    bool test_columnPermute(comp_r_mat* A, int col1, int col2);
    
    bool test_changeElement( comp_r_mat* A , int rowInd , int colInd , double newValue );
    
    bool test_scalarMultiple( comp_r_mat* A , double scale );
    
    bool test_copyMatrix( comp_r_mat* C , comp_r_mat* A );
    
    bool test_decomposeMatrix( vector<double>* DS , comp_r_mat* LUS , comp_r_mat* AS );
    
    bool test_matrixProduct ( vector<double>* result , comp_r_mat* A , vector<double>* vec);
    
    bool test_jacobiSolver( vector<double>* X , vector<double>* DS , comp_r_mat* LUS , vector<double>* B );
    
    bool test_calculateNorm( double& norm , vector<double>* v , vector<double>* Ax );
}


#endif /* test_comp_func_hpp */
