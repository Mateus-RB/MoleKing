//   MoleKing //
//
//   File:        [Matrix.cpp]
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

#include "Matrix.hpp"
#include <string>

Matrix::Matrix(vector < vector <double> > matrix){
    this->matrix = matrix;
};

Matrix::Matrix(){};

void Matrix::setMatrix(vector < vector <double> > matrix){
    this->matrix = matrix;
}

Matrix::~Matrix(){
    this->matrix.clear();
    this->matrix.resize(0);
};

Matrix::Matrix(int i, int j){
    this->matrix = vector < vector < double > > (i, vector < double > (j));
};

void Matrix::replace(int i, int j, double newValue){
    this->matrix.at(i).at(j) = newValue;
};

Matrix Matrix::sum(Matrix matrix_B){
    vector <vector <double> > matrixB = matrix_B.toVector();
    vector < vector < double > > result(this->matrix.size(), vector <double> (this->matrix[0].size()));
    if (this->matrix.size() == matrixB.size() && this->matrix[0].size() == matrixB[0].size()){
        for (int i = 0; i < (int) this->matrix.size(); i++){
            for(int j = 0; j < (int) this->matrix[0].size(); j++){
                double newValue = this->matrix[i][j] + matrixB[i][j];
                result.at(i).at(j) = newValue;
            };
        };
    } else {
        exit (EXIT_FAILURE);
    };
    Matrix m_result = Matrix(result);
    return m_result;
};


Matrix Matrix::multiplication(Matrix matrix_B){
    vector <vector <double> > matrixB = matrix_B.toVector();
    vector <double> in(matrixB[0].size(), 0);
    vector < vector < double > > result(this->matrix.size(), in);
    if (this->matrix[0].size() == matrixB.size()){
        for (int i = 0; i < (int) this->matrix.size(); i++){
            for (int j = 0; j < (int) matrixB[0].size(); j++){
                double newValue = 0;
                for(int k = 0; k < (int) this->matrix[0].size(); k++){
                    newValue = newValue + this->matrix[i][k] * matrixB[k][j];
                };
                result.at(i).at(j) = newValue;
            };
        };
    } else {
        exit (EXIT_FAILURE);
    };
    Matrix m_result = Matrix(result);
    return m_result;
};

Matrix Matrix::multiplication(double scalar){
    vector < vector < double > > result(this->matrix.size(), vector <double> (this->matrix[0].size()));
    for (int i = 0; i < (int) this->matrix.size(); i++){
        for (int j = 0; j < (int) this->matrix[0].size(); j++){
            result.at(i).at(j) = this->matrix.at(i).at(j) * scalar;
        };
    };
    Matrix m_result = Matrix(result);
    return m_result;
};

vector <long> Matrix::getDimensions(){
    long lines = this->matrix.size();
    long columns = this->matrix[0].size();
    return vector <long> {lines, columns};
};

double Matrix::determinant(){
    vector <long> dimensions = this->getDimensions();
    if(dimensions[0] != dimensions[1]){
        exit(EXIT_FAILURE);
    };
    long n = dimensions[0];
    double D;
    return D = det(this->matrix, n); 
};

double Matrix::det(vector< vector <double> > mat, long n){
    double D = 0;
    if (n == 1){
        return mat.at(0).at(0); 
    }
    vector < vector <double> > temp;
  
    int sign = 1;

    for (int f = 0; f < n; f++){ 
        vector < vector <double> > temp = getCofactor(mat, 0, f, n);
        D += sign * mat.at(0).at(f) * this->det(temp, n - 1); 
        sign = -sign; 
    };
    return D; 
};

vector < vector <double> >  Matrix::getCofactor(vector< vector <double> > mat, int p, int q, long n){
    vector <double> line(n-1);
    vector < vector <double> > temp(n-1, line);
    int i = 0, j = 0; 
    for (int row = 0; row < n; row++) { 
        for (int col = 0; col < n; col++){ 
            if (row != p && col != q) { 
                temp.at(i).at(j++) = mat.at(row).at(col);
                if (j == n - 1) { 
                    j = 0; 
                    i++; 
                };
            }; 
        }; 
    };
    return temp;
};

double Matrix::element(long i = 1, long j = 1){
    return matrix[i-1][j-1];
};

void Matrix::print(){
    for (int i = 0; i < (int) this->matrix.size(); i++){
        for (int j = 0; j < (int) this->matrix[0].size(); j++){
            cout << this->matrix[i][j] << " ";
        }
        cout << endl;
    }
};

vector <double> Matrix::getLine(int i){
    return this->matrix[i];
};

string Matrix::toStr(){
    string temp;
    for (int i = 0; i < (int) this->matrix.size(); i++){
        for (int j = 0; j < (int) this->matrix[0].size(); j++){
            temp = temp + std::to_string(this->matrix[i][j]) + " ";
        };
        temp = temp + "\n";
    };
    return temp;
};

vector < vector <double> > Matrix::toVector(){
    return this->matrix;
}
