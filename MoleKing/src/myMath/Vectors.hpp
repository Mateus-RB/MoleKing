//   MoleKing //
//
//   File:        [Vectors.hpp]
//
//   Author(s):   ['LEEDMOL Reserch Group']
//   Site(s):     ['https://www.researchgate.net/lab/LEEDMOL-Heibbe-Cristhians-Research-Group-Heibbe-Cristhian-B-De-Oliveira']
//   Email(s):    ['heibbe@ufg.br']
//   Credits:     ['Copyright © 2023 LEEDMOL. All rights reserved.']
//   Date:        ['17.01.2023']
//   Version:     ['1.4.2']
//   Status:      ['Development']
//   Language:    ['C++','Python']
//   Description: ['A python module written in C++ for theoretical chemestry']


#ifndef Vectors_hpp
#define Vectors_hpp

#include <stdio.h>
#include <vector>
#include <string>
#include <math.h>
#include <iostream>

using namespace std;

class Vector3D{
    private:
    double x_a, x_b, y_a, y_b, z_a, z_b, s_i, s_j, s_k;

    public:
    Vector3D(vector<double> pointA, vector<double> pointB = {0.0, 0.0, 0.0});
    Vector3D();
    void setVector(vector<double> pointA, vector<double> pointB = {0.0, 0.0, 0.0});
    ~Vector3D();
    double magnitude();
    vector <double> getVector();
    void show();
    Vector3D normalize();
    Vector3D conjugate();
    Vector3D operator/ (double mag);
    Vector3D operator* (double mag);
    Vector3D operator+ (Vector3D vectorB);
    Vector3D operator- (Vector3D vectorB);
    Vector3D crossProduct(Vector3D vectorB);
    double dotProduct(Vector3D vectorB);
    double angle(Vector3D vectorB, char unit = 'd');
    double axisValue(char unitVector);
    string toStr();
};

class Quaternion{
    private:
    double u, s_i, s_j, s_k;

    public:
    Quaternion(double u, vector <double> vectorA, vector <double> vectorB);
    ~Quaternion();
    double magnitude();
    vector <double> getQuaternion();
    void show();
};

#endif /* Vectors_hpp */
