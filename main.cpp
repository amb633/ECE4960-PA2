#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>

#include "compressed_func.hpp"

using namespace std;
using namespace compressed;

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

    const int rank = row_ptr.size() - 1;

    //1: create row-compressed matrix
    comp_r_mat AC;
    AC.value = values;
    AC.col_id = cols;
    AC.row_p = row_ptr;

    //2: decompose matrix to diagonals and LU form
    comp_r_mat LUC;
    vector<double> DC;
    decomposeMatrix( &DC , &LUC , &AC );

    //3: create first B vector and solution vector
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
    calculateNorm( normB , &B , &zeros );

    //5: first solver loop
    cout << "---------- first solver loop ---------- " << endl;
    while( abs(normCurrent - normPrev) > 1e-10 ){
        normPrev = normCurrent;
        jacobiSolver( &X , &DC , &LUC , &B );
        productAx( &matProd , &AC , &X );
        calculateNorm( normCurrent , &B , &matProd );
        cout << counter << " : ";
        /*for ( int i = 0 ; i < 5 ; i++ ){
            cout << X[i] << "   " ;
        }*/
        cout << "  residual norm = " << normCurrent/normB << endl;
        counter++;
    }
    cout << endl;

    //6: second solver loop
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
    calculateNorm( normB , &B , &zeros );

    cout << "---------- second solver loop ---------- " << endl;
    while( abs(normCurrent - normPrev) > 1e-10 ){
        normPrev = normCurrent;
        jacobiSolver( &X , &DC , &LUC , &B );
        productAx( &matProd , &AC , &X );
        calculateNorm( normCurrent , &B , &matProd );
        cout << counter << " : ";
        /*for ( int i = 0 ; i < 5 ; i++ ){
            cout << X[i] << "   " ;
        }*/
        cout << "  residual norm = " << normCurrent/normB << endl;
        counter++;
    }
    cout <<endl;

    //6: third solver loop
    for ( int i = 0 ; i < rank ; i++ ){
        X[i] = 1.0/DC[i];
        matProd[i] = 0.0;
        B[i] = 1.0;
    }

    normPrev = 2;
    normCurrent = 1;
    normB = 0;
    counter = 0;
    calculateNorm( normB , &B , &zeros );
    
    cout << "---------- third solver loop ---------- " << endl;
    while( abs(normCurrent - normPrev) > 1e-10 ){
        normPrev = normCurrent;
        jacobiSolver( &X , &DC , &LUC , &B );
        productAx( &matProd , &AC , &X );
        calculateNorm( normCurrent , &B , &matProd );
        cout << counter << " : ";
        /*for ( int i = 0 ; i < 5 ; i++ ){
            cout << X[i] << "   " ;
        }*/
        cout << "  residual norm = " << normCurrent/normB << endl;
        counter++;
    }
    cout << endl;

    return 0;
}