//
//  Matrix.hh
//  MoleKing_util
//
//  Created by Thiago Lopes, Sandro Brito and Mateus Barbosa on 14/01/20.
//  Copyright © 2020 LMSC. All rights reserved.
//

#ifndef Matrix_hpp
#define Matrix_hpp
#include <math.h>
#include <vector>


#include <iostream>

using namespace std;

class Matrix{
    private:
    vector < vector <double> > matrix;
    vector < vector <double> > getCofactor(vector< vector <double> > mat, int p, int q, long n);
    double det(vector< vector <double> > mat, long n);

    public:
    Matrix(vector < vector <double> > matrix);
    Matrix();
    Matrix(int i, int j);
    ~Matrix();
    void setMatrix(vector < vector <double> > matrix);
    Matrix sum(Matrix matrixB);
    Matrix multiplication(double scalar);
    vector < vector < double > > toVector();
    Matrix multiplication(Matrix matrixB);
    double determinant();
    void replace(int i, int j, double newValue);
    vector<double> getLine(int i);
    vector <long> getDimensions();
    double element(long i, long j);
    void print();
    string toStr();
};

#endif /* Matrix_hpp */
