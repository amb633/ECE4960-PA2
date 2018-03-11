#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <chrono>

#include "compressed_func.hpp"
#include "full_func.hpp"
#include "test_comp_func.hpp"
#include "test_full_func.hpp"
#include "wilkinson_test_func.hpp"
#include "large_matrix_test.hpp"

using namespace std;

int main(int argc, char const *argv[])
{
    cout << endl;
    cout << "-------------------- testing full matrix functions -------------------- " << endl;   
    test_full::run_full_matrix_tests();
    cout << endl << endl;
    cout << "-------------------- testing compressed matrix functions -------------------- " << endl;
    test_compressed::call_comp_tests();
    cout << endl << endl;
    cout << "-------------------- wilkinson tests across full and compressed matrices -------------------- " << endl;
    wilkinson_test::run_wilkinson_tests();
    cout << endl << endl;
    cout << "-------------------- large matrix tests -------------------- " << endl;
    large_matrix_test::run_large_matrix_tests();
    cout << endl << endl;
    
    auto start_creation = chrono::system_clock::now();
    
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

    const int rank = row_ptr.size() - 1;

    //1: create row-compressed matrix
    compressed::comp_r_mat AC;
    AC.value = values;
    AC.col_id = cols;
    AC.row_p = row_ptr;

    auto end_creation = chrono::system_clock::now();
    cout << endl;

    //2: decompose matrix to diagonals and LU form
    compressed::comp_r_mat LUC;

    vector<double> DC;
    compressed::decomposeMatrix( &DC , &LUC , &AC );

    //3: create first B vector and solution vector
    auto first_pre_start = chrono::system_clock::now();

    vector<double> B , X , matProd , zeros ;
    for ( int i = 0 ; i < rank ; i++ ){
        B.push_back(0.0);
        X.push_back(0.0);
        matProd.push_back(0.0);
        zeros.push_back(0.0);
    }
    B[0] = 1.0;
    X[0] = (1.0/DC[0]);

    //4: initialize variables for loop
    double normPrev = 2;
    double normCurrent = 1;
    double normB = 0;
    int counter = 0;
    compressed::calculateNorm( normB , &B , &zeros );

    //5: first solver loop
    auto first_start = chrono::system_clock::now();
    cout << "-------------------- first solver loop -------------------- " << endl;
    while( abs(normCurrent - normPrev) > 1e-7 ){
        normPrev = normCurrent;
        compressed::jacobiSolver( &X , &DC , &LUC , &B );
        compressed::productAx( &matProd , &AC , &X );
        compressed::calculateNorm( normCurrent , &B , &matProd );
        cout << counter << " : ";
        /*for ( int i = 0 ; i < 5 ; i++ ){
            cout << X[i] << "   " ;
        }*/
        cout << "  residual norm = " << normCurrent/normB << endl;
        counter++;
    }
    cout << endl;
    auto first_end = chrono::system_clock::now();

    //6: second solver loop
    auto second_pre_start = chrono::system_clock::now();
    for ( int i = 0 ; i < rank ; i++ ){
        X[i] = 0.0;
        matProd[i] = 0.0;
    }
    B[0] = 0.0;
    B[4] = 1.0;
    X[4] = (1.0/DC[4]);

    normPrev = 2;
    normCurrent = 1;
    normB = 0;
    counter = 0;
    compressed::calculateNorm( normB , &B , &zeros );

    cout << "-------------------- second solver loop -------------------- " << endl;
    auto second_start = chrono::system_clock::now();
    while( abs(normCurrent - normPrev) > 1e-7 ){
        normPrev = normCurrent;
        compressed::jacobiSolver( &X , &DC , &LUC , &B );
        compressed::productAx( &matProd , &AC , &X );
        compressed::calculateNorm( normCurrent , &B , &matProd );
        cout << counter << " : ";
        /*for ( int i = 0 ; i < 5 ; i++ ){
            cout << X[i] << "   " ;
        }*/
        cout << "  residual norm = " << normCurrent/normB << endl;
        counter++;
    }
    cout <<endl;
    auto second_end = chrono::system_clock::now();

    //6: third solver loop
    auto third_pre_start = chrono::system_clock::now();
    for ( int i = 0 ; i < rank ; i++ ){
        X[i] = 1.0/DC[i];
        matProd[i] = 0.0;
        B[i] = 1.0;
    }

    normPrev = 2;
    normCurrent = 1;
    normB = 0;
    counter = 0;
    compressed::calculateNorm( normB , &B , &zeros );
    
    cout << "-------------------- third solver loop -------------------- " << endl;
    auto third_start = chrono::system_clock::now();
    while( abs(normCurrent - normPrev) > 1e-7 ){
        normPrev = normCurrent;
        compressed::jacobiSolver( &X , &DC , &LUC , &B );
        compressed::productAx( &matProd , &AC , &X );
        compressed::calculateNorm( normCurrent , &B , &matProd );
        cout << counter << " : ";
        /*for ( int i = 0 ; i < 5 ; i++ ){
            cout << X[i] << "   " ;
        }*/
        cout << "  residual norm = " << normCurrent/normB << endl;
        counter++;
    }
    cout << endl;
    auto third_end = chrono::system_clock::now();

    chrono::duration<double> elapsed_creation = end_creation - start_creation;
    chrono::duration<double> elapsed_first_pre_start = first_start - first_pre_start;
    chrono::duration<double> elapsed_first_solver = first_end - first_start;
    chrono::duration<double> elapsed_second_pre_start = second_start - second_pre_start;
    chrono::duration<double> elapsed_second_solver = second_end - second_start;
    chrono::duration<double> elapsed_third_pre_start = third_start - third_pre_start;
    chrono::duration<double> elapsed_third_solver = third_end - third_start;

    cout << "Time taken for reading file and creating large matrix = " << elapsed_creation.count() << endl;
    cout << "Time taken for initializing variables for first solver = " << elapsed_first_pre_start.count() << endl;
    cout << "Time taken for first solver loop = " << elapsed_first_solver.count() << endl;
    cout << "Time taken for initializing variables for second solver= " << elapsed_second_pre_start.count() << endl;
    cout << "Time take for second solver loop = " << elapsed_second_solver.count() << endl;
    cout << "Time taken for initializing variables for third solver = " << elapsed_third_pre_start.count() << endl;
    cout << "Time taken for third solver loop = " << elapsed_third_solver.count() << endl;

    return 0;
}
