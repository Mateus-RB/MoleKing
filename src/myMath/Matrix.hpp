//   MoleKing //
//
//   File:        [Matrix.hpp]
//
//   Author(s):   ['LEEDMOL Reserch Group']
//   Site(s):     ['https://www.researchgate.net/lab/LEEDMOL-Heibbe-Cristhians-Research-Group-Heibbe-Cristhian-B-De-Oliveira']
//   Email(s):    ['heibbe@ufg.br']
//   Credits:     ['Copyright © 2023 LEEDMOL. All rights reserved.']
//   Date:        ['17.01.2023']
//   Version:     ['1.5.1']
//   Status:      ['Development']
//   Language:    ['C++','Python']
//   Description: ['A python module written in C++ for theoretical chemestry']


#ifndef Matrix_hpp
#define Matrix_hpp
#include <cmath>
#include <vector>
//#include <../../include/Eigen/Core>


#include <iostream>

using namespace std;

class Matrix{
    private:
    vector < vector <double> > matrix;
    vector < vector <double> > getCofactor(vector< vector <double> > mat, int p, int q, long n);
    double det(vector< vector <double> > mat, long n);

    public:

    //double eigenvalue();

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
