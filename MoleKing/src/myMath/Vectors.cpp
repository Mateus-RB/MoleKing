//   MoleKing //
//
//   File:        [Vectors.cpp]
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

#include "Vectors.hpp"

// Vector3D class ///

Vector3D::Vector3D(vector<double> pointA, vector<double> pointB){
    this->x_a = pointA[0];
    this->x_b = pointB[0];
    this->s_i = this->x_a - this->x_b;
    this->y_a = pointA[1];
    this->y_b = pointB[1];
    this->s_j = this->y_a - this->y_b;
    this->z_a = pointA[2];
    this->z_b = pointB[2];
    this->s_k = this->z_a - this->z_b;
};

Vector3D::Vector3D(){};

string Vector3D::toStr(){
    string temp;
    temp = to_string(this->s_i) + "i";
    if(this->s_j >= 0){
        temp = temp + " +" + to_string(this->s_j) + "j";
    } else {
        temp = temp + " -" + to_string(this->s_j) + "j";
    };
    if(this->s_k >= 0){
        temp = temp + " +" + to_string(this->s_k) + "k";
    } else {
        temp = temp + " -" + to_string(this->s_k) + "k";
    };
    return temp;
}

void Vector3D::setVector(vector<double> pointA, vector<double> pointB){
    this->x_a = pointA[0];
    this->x_b = pointB[0];
    this->s_i = this->x_a - this->x_b;
    this->y_a = pointA[1];
    this->y_b = pointB[1];
    this->s_j = this->y_a - this->y_b;
    this->z_a = pointA[2];
    this->z_b = pointB[2];
    this->s_k = this->z_a - this->z_b;
};

Vector3D::~Vector3D(){
    x_a = 0.0;
    x_b = 0.0;
    y_a = 0.0;
    y_b = 0.0;
    z_a = 0.0;
    z_b = 0.0;
    s_i = 0.0;
    s_j = 0.0;
    s_k = 0.0;
};

double Vector3D::magnitude(){
    double norm = pow(this->s_i, 2) + pow(this->s_j, 2) + pow(this->s_k, 2);
    return sqrt(norm);
};

vector <double> Vector3D::getVector(){
    return vector <double> {this->s_i, this->s_j, this->s_k};
};

Vector3D Vector3D::normalize(){
    Vector3D v = *this;
    double mag = 1/this->magnitude();
    return Vector3D(vector <double> {this->s_i * mag, this->s_j * mag, this->s_k * mag});
};

Vector3D Vector3D::conjugate(){
    Vector3D v = *this;
    double mag = -1.0;
    return Vector3D(vector <double> {this->s_i * mag, this->s_j * mag, this->s_k * mag});
};

Vector3D Vector3D::operator/ (double mag){
    Vector3D v = *this;
    return Vector3D(vector <double> {this->s_i / mag, this->s_j / mag, this->s_k / mag});
};

Vector3D Vector3D::operator* (double mag){
    Vector3D v = *this;
    return Vector3D(vector <double> {this->s_i * mag, this->s_j * mag, this->s_k * mag});
};

Vector3D Vector3D::crossProduct(Vector3D vectorB){
    Vector3D v = *this;
    vector < double > b = vectorB.getVector();
    return Vector3D(vector <double> {(this->s_j * b[2]) - (this->s_k * b[1]), (this->s_k * b[0]) - (this->s_i * b[2]), (this->s_i *b[1]) - (b[0] * this->s_j)});
};

double Vector3D::dotProduct(Vector3D vectorB){
    Vector3D v = *this;
    vector < double > b = vectorB.getVector();
    double r = double (this->s_i * b[0] + this->s_j * b[1] + this->s_k * b[2]);
    return r;
};

Vector3D Vector3D::operator+ (Vector3D vectorB){
    Vector3D v = *this;
    vector < double > b = vectorB.getVector();
    return Vector3D(vector <double> {this->s_i + b[0], this->s_j + b[1], this->s_k + b[2]});
};

Vector3D Vector3D::operator-(Vector3D vectorB){
    Vector3D v = *this;
    return this->operator+(vectorB.conjugate());
};

double Vector3D::angle(Vector3D vectorB, char unit){
    double tetha;

    tetha = acos( this->dotProduct(vectorB) / ( vectorB.magnitude() * this->magnitude() ) );
    if (unit == 'd'){
        return (tetha * 180) / M_PI;
    } else {
        return tetha;
    };
        
};

double Vector3D::axisValue(char unitVector){
    if (unitVector == 'i' || unitVector == 'x'){
        return this->s_i;
    } else if (unitVector == 'j' || unitVector == 'y'){
        return this->s_j;
    } else if (unitVector == 'k' || unitVector == 'z'){
        return this->s_k;
    } else {
        return 0.0;
    };
};

// Quaternion class //

Quaternion::Quaternion(double u, vector <double> vectorA, vector <double> vectorB = {0.0, 0.0, 0.0}){
    //Quaternion q = *this;
    this->u = u;
    this->s_i = vectorA[0] - vectorB[0];
    this->s_j = vectorA[0] - vectorB[0];
    this->s_k = vectorA[0] - vectorB[0];
};

Quaternion::~Quaternion(){
    u = 0.0;
    s_i = 0.0;
    s_j = 0.0;
    s_k = 0.0;
};

double Quaternion::magnitude(){
    double norm = pow(this->u, 2) + pow(this->s_i, 2) + pow(this->s_j, 2) + pow(this->s_k, 2);
    return sqrt(norm);
};

void Quaternion::show(){
    cout << "q = " << this->u << " + " << this->s_i << "i +" << this->s_j << "j +" << this->s_k << "k" << endl;
}

vector <double> Quaternion::getQuaternion(){
    return vector <double> {this->u, this->s_i, this->s_j, this->s_k};
};
