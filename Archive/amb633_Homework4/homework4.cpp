//
//  main.cpp
//  Week4
//
//  Created by Ariana Bruno on 2/13/18.
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
    comp_r_mat mat_a = construct_compressed_matrix(input_vector);
    comp_r_mat mat_b = construct_compressed_matrix(single_vector);
    print_comp_r_mat( &mat_a );
    comp_r_mat reorder_a = mat_a;
    comp_r_mat reorder_b = mat_b;
    
    // assumes a sqaure matrix
    for( int k = 0; k< mat_a.row_p.size()-1; k++ ){
        reorderMat(&reorder_a, &reorder_a, &reorder_b, k, k);
        print_comp_r_mat(&reorder_a);
        print_comp_r_mat(&reorder_b);
    }

// Information below is from previous assignment not relevant at the moment
    
//    comp_r_mat mat_a1 = mat_a;
//    rowPermute(&mat_a1, 0, 2);
//    cout << "After the row permutaiton of 1 and 3" << endl;
//    print_comp_r_mat( &mat_a1 );
//
//    rowPermute(&mat_a1, 0, 4);
//    cout << "After the row permutaiton of 1 and 5" << endl;
//    print_comp_r_mat( &mat_a1 );
//
//    comp_r_mat mat_a2 = mat_a;
//    rowScale(&mat_a2, 0, 3, 3.0);
//    cout << "After the row scaling of 3*row[1] + row[4]" << endl;
//    print_comp_r_mat( &mat_a2 );
//
//
//    cout << "Testing the A*x + b: " << endl;
//    vector<double> x_test = {5,4,3,2,1};
//    vector<double> b_test;
//    int check_product = productAx(&mat_a, &x_test, &b_test);
//    cout << check_product << endl;
//    if (check_product == 0) {
//        cout << "This is the resulting product of A*x = b, b= { ";
//        for (int i = 0; i<b_test.size(); i++){
//            cout << b_test[i] << " ";
//        }
//        cout << "} " << endl;
//    }
//
//
//    ifstream memplus("/Users/arianabruno/Desktop/ECE4960/memplus.mtx");
//    for(int i = 0; i<1; i++){
//        memplus.ignore(2048, '\n');
//    }
//    int colRank, rowRank, maxEntry;
//    memplus >> colRank >> rowRank >> maxEntry;
//    int i, j;
//    double v;
//    vector<int> I, J;
//    vector<double> Val;
//    while( memplus >> i >> j >> v ){
//        I.push_back(i);
//        J.push_back(j);
//        Val.push_back(v);
//    };
//    comp_r_mat mat_memplus = construct_compressed_matrix( &I, &J, &Val, rowRank, colRank);
//
//    comp_r_mat mem1 = mat_memplus;
//    rowPermute(&mem1, 1, 3);
//    rowPermute(&mem1, 1, 5);
//    rowPermute(&mem1, 10, 3000);
//    rowPermute(&mem1, 5000, 10000);
//    rowPermute(&mem1, 6, 15000);
//
//    comp_r_mat mem2 = mat_memplus;
//    rowScale(&mem2, 2, 4, 3.0);
//    rowPermute(&mem2, 2, 5);
//    rowScale(&mem2, 5, 4, -3.0);
//    vector<double> x_test2(rowRank+1, 1.0);
//    vector<double> b_test2;
//    int check_product2 = productAx(&mat_memplus, &x_test2, &b_test2);
//    cout << "This is the resulting to check product of memplus * vec_x(of all 1's) = b" << endl;
//    cout << check_product2 << endl;
//    cout<< boolalpha;
//    cout << "This is the check sum results: " << check_sum(&mat_memplus, &b_test2) << endl;

    return 0;
}


