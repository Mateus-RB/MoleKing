//   MoleKing //
//
//   File:        [Geometry.cpp]
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


#include "Geometry.hpp"

// Point class //

Point::Point(){};

Point::Point(double coord1, double coord2, double coord3, char typeCoord){
    if (typeCoord == 'c'){
        this->x = coord1;
        this->y = coord2;
        this->z = coord3;
        vector<double> temp = SphericalCoords(this->x, this->y, this->z, 'c').toSpherical();
        this->radius = temp[0];
        this->tetha = temp[1];
        this->phi = temp[2];
    } else {
        this->radius = coord1;
        this->tetha = coord2;
        this->phi = coord3;
        vector<double> temp = SphericalCoords(this->radius, this->tetha, this->phi, 's').toCartesian();
        this->x = temp[0];
        this->y = temp[1];
        this->z = temp[2];
    };
};

void Point::setPoint(double coord1, double coord2, double coord3, char typeCoord){
    if (typeCoord == 'c'){
        this->x = coord1;
        this->y = coord2;
        this->z = coord3;
        vector<double> temp = SphericalCoords(this->x, this->y, this->z, 'c').toSpherical();
        this->radius = temp[0];
        this->tetha = temp[1];
        this->phi = temp[2];
    } else {
        this->radius = coord1;
        this->tetha = coord2;
        this->phi = coord3;
        vector<double> temp = SphericalCoords(this->radius, this->tetha, this->phi, 's').toCartesian();
        this->x = temp[0];
        this->y = temp[1];
        this->z = temp[2];
    };
};

Point::~Point(){
    this->x = 0.0;
    this->y = 0.0;
    this->z = 0.0;
    this->radius = 0.0;
    this->tetha = 0.0;
    this->phi =  0.0;
};

bool Point::operator==(Point point2){
    bool value;
    if (this->x == point2.getCoords('c')[0]){
        if (this->y == point2.getCoords('c')[1]){
            if (this->z == point2.getCoords('c')[2]){
                value = 1;
            } else {
                value = 0;
            };
        } else {
            value = 0;
        };
    } else {
        value = 0;
    };
    return value;
};

void Point::setCoords(vector <double> newValues, char typeCoord = 'c' /* 'c' for cartesian coordinates, 's' for spherical*/){
    if(typeCoord == 'c'){
        this->x = newValues[0];
        this->y = newValues[1];
        this->z = newValues[2];
        vector<double> temp = SphericalCoords(this->x, this->y, this->z, 'c').toSpherical();
        this->radius = temp[0];
        this->tetha = temp[1];
        this->phi = temp[2];
    } else {
        this->radius = newValues[0];
        this->tetha = newValues[1];
        this->phi = newValues[2];
        vector<double> temp = SphericalCoords(this->radius, this->tetha, this->phi, 's').toCartesian();
        this->x = temp[0];
        this->y = temp[1];
        this->z = temp[2];
    };
};

void Point::setCoord(char coordName, double newValue){
    if(coordName == 'x'){
        this->x = newValue;
    } else if (coordName == 'y'){
        this->y = newValue;
    } else if(coordName == 'z'){
        this->z = newValue;
    } else if(coordName == 'r'){
        this->radius = newValue;
    } else if(coordName == 't'){
        this->tetha = newValue;
    } else if(coordName == 'p'){
        this->phi = newValue;
    } else {
        exit (EXIT_FAILURE);
    };
    if (coordName == 'x' || coordName == 'y' || coordName == 'z'){
        this->setCoords(vector <double> {this->x, this->y, this->z}, 'c');
    } else {
        this->setCoords(vector <double> {this->radius, this->tetha, this->phi}, 's');
    };
};

vector<double> Point::getCoords(char typeCoord = 'c'){
    if (typeCoord == 'c'){
        return vector <double> {this->x, this->y, this->z};
    } else {
        return vector <double> {this->radius, this->tetha, this->phi};
    };
};

void Point::translation(Vector3D translationVector){
    vector < vector < double> > posMatrix= { {this->x}, {this->y}, {this->z}, {1.0} };
    Matrix transMAtrix = Matrix( { {1.0, 0.0, 0.0, translationVector.axisValue('i')},
                                   {0.0, 1.0, 0.0, translationVector.axisValue('j')},
                                   {0.0, 0.0, 1.0, translationVector.axisValue('k')},
                                   {0.0, 0.0, 0.0, 1.0} } );
    Matrix newPos = transMAtrix.multiplication(posMatrix);
    this->x = newPos.element(1, 1);
    this->y = newPos.element(2, 1);
    this->z = newPos.element(3, 1);
    vector <double> newPosSpherical = SphericalCoords(this->x, this->y, this->z, 'c').toSpherical();
    this->radius = newPosSpherical[0];
    this->tetha = newPosSpherical[1];
    this->phi = newPosSpherical[2];
};


void Point::rotationVector(double angle, Vector3D unitVector){
    angle = M_PI * angle / 180;
    double x_u = unitVector.axisValue('i');
    double y_u = unitVector.axisValue('j');
    double z_u = unitVector.axisValue('k');
    vector < vector < double> > posMatrix= { {this->x}, {this->y}, {this->z}, {1.0} };

    Matrix rotMatrix = Matrix( { { cos(angle) + pow(x_u, 2)*(1 - cos(angle)), y_u * x_u * (1 - cos(angle)) - z_u * (sin(angle)), z_u * x_u * (1 - cos(angle)) + y_u * (sin(angle)), 0.0},
                                 { x_u * y_u * (1 - cos(angle)) + z_u * (sin(angle)), cos(angle) + (1 - cos(angle)) * pow(y_u, 2), z_u * y_u * (1 - cos(angle)) - x_u * (sin(angle)), 0.0},
                                 { x_u * z_u * (1 - cos(angle)) - y_u * sin(angle), y_u * z_u * (1 - cos(angle)) + x_u * sin(angle), cos(angle) + (1 - cos(angle)) * pow(z_u, 2), 0.0},
                                 {0.0, 0.0, 0.0, 1.0} } );

    Matrix newPos = rotMatrix.multiplication(posMatrix);
    this->x = newPos.element(1, 1);
    this->y = newPos.element(2, 1);
    this->z = newPos.element(3, 1);
    vector <double> newPosSpherical = SphericalCoords(this->x, this->y, this->z, 'c').toSpherical();
    this->radius = newPosSpherical[0];
    this->tetha = newPosSpherical[1];
    this->phi = newPosSpherical[2];
}

string Point::toStr(char spaceType){
    string temp = "Coords in ";
    if (spaceType == 'c'){
        temp = temp + "Cartesian Space (x, y, z): ";
        temp = temp + "(" + to_string(this->x) + ", " + to_string(this->y) + ", " + to_string(this->z) + ")";
    } else {
        temp = temp + "Spherical Space: (Radius, Polar Angle, Azimuthal Angle)";
        temp = temp + "(" + to_string(this->radius) + ", " + to_string(this->tetha) + ", " + to_string(this->phi) + ")";
    };
    return temp;
};

// SphericalCoords class //

SphericalCoords::SphericalCoords(double coord1/*x or radius*/, double coord2/*y or tetha*/, double coord3/*z or phi*/, char spaceType = 'c'/* 'c' for cartesian ou 's' for spherical*/){
    if (spaceType == 'c'){
        this->x = coord1;
        this->y = coord2;
        this->z = coord3;
    } else {
        this->radius = coord1;
        this->tetha =  coord2;
        this->phi = coord3;
    };
};

SphericalCoords::~SphericalCoords(){
    this->x = 0.0;
    this->y = 0.0;
    this->z = 0.0;
    this->radius = 0.0;
    this->tetha = 0.0;
    this->phi = 0.0;
};

vector <double> SphericalCoords::toCartesian(){
    SphericalCoords s = *this;
    s.x = s.radius * sin(M_PI * s.tetha / 180) * cos(M_PI * s.phi / 180);
    s.y = s.radius * sin(M_PI * s.tetha / 180) * sin(M_PI * s.phi / 180);
    s.z = s.radius * cos(M_PI * s.tetha / 180);
    return vector <double> {s.x, s.y, s.z};
};

vector <double> SphericalCoords::toSpherical(){
    SphericalCoords c = *this;
    c.radius = sqrt(pow(c.x, 2) + pow(c.y, 2) + pow(c.z, 2));
    if (c.radius == 0){
        return vector <double> {0, 0, 0};
    }
    c.tetha = acos(c.z/c.radius) * 180 / M_PI;
    double xy = sqrt(pow(c.x, 2) + pow(c.y,2));
    if (xy == 0){
        c.phi =0;
    } else {
        c.phi = acos(c.x/xy) * 180 / M_PI;
    };
    return vector <double> {c.radius, c.tetha, c.phi};
};

// StraightSegment class //

StraightSegment::StraightSegment(Point a, Point b){
    this->a = a;
    this->b = b;
    this->calcAbs();
};

double StraightSegment::getValue(){
    return this->absValue;
};

vector <Point> StraightSegment::getPoints(){
    return vector <Point> {this->a, this->b};
};

void StraightSegment::calcAbs(){
    this->absValue = Vector3D(this->a.getCoords(), this->b.getCoords()).magnitude();
};

void StraightSegment::stretchNcontract(double increment, char freezePoint){
    double factor = (this->absValue + increment)/this->absValue;
    vector<double> oldA = this->a.getCoords();
    vector<double> oldB = this->b.getCoords();
    if (freezePoint == 'a'){
        Vector3D vec = Vector3D(oldB, oldA);
        vec = vec * factor;
        vector <double> vTrans = vec.getVector();
        this->a = Point(oldA[0], oldA[1],oldA[2]);
        this->b = Point(vTrans[0]+oldA[0], vTrans[1]+oldA[1], vTrans[2]+oldA[2]);
        this->calcAbs();
    } else if (freezePoint == 'b'){
        Vector3D vec = Vector3D(oldA, oldB);
        vec = vec * factor;
        vector <double> vTrans = vec.getVector();
        this->a = Point(vTrans[0]+oldB[0], vTrans[1]+oldB[1], vTrans[2]+oldB[2]);
        this->b = Point(oldB[0], oldB[1], oldB[2]);
        this->calcAbs();
    } else{
        cout << "stretchNcontract method of StraightSegment." << endl;
        exit(0);
    }
};

// Angle class //

Angle::Angle(Point a, Point b, Point c){
    this->a = a;
    this->b = b;
    this->c = c;
    this->calcAbs();
};

void Angle::calcAbs(){
    Vector3D r1 = Vector3D(this->a.getCoords(), this->b.getCoords());
    Vector3D r2 = Vector3D(this->c.getCoords(), this->b.getCoords());
    this->absValue = r1.angle(r2);
};

double Angle::getValue(){
    return this->absValue;
};

vector <Point> Angle::getPoints(){
    return vector <Point> {this->a, this->b, this->c};
};

void Angle::increaseNdecrease(double increment, char freezePoint){
    Vector3D r1 = Vector3D(this->a.getCoords(), this->b.getCoords());
    Vector3D r2 = Vector3D(this->c.getCoords(), this->b.getCoords());
    Vector3D norm = r1.crossProduct(r2);
    if (freezePoint == 'a'){
        this->a.rotationVector(this->absValue, norm);
    } else {
        this->c.rotationVector(this->absValue, norm);
    }
    this->calcAbs();
};

// Torsion class //

Torsion::Torsion(Point a, Point b, Point c, Point d){
    this->a = a;
    this->b = b;
    this->c = c;
    this->d = d;
    this->calcAbs();
};

void Torsion::calcAbs(){
    Vector3D r1 = Vector3D(this->b.getCoords(), this->a.getCoords());
    Vector3D r2 = Vector3D(this->b.getCoords(), this->c.getCoords());
    Vector3D r3 = Vector3D(this->c.getCoords(), this->d.getCoords());
    Vector3D semi_normal1 = r1.crossProduct(r2) / sin(r1.angle(r2, 'r'));
    Vector3D semi_normal2 = r3.crossProduct(r2) / sin(r3.angle(r2, 'r'));
    double angleD = semi_normal1.angle(semi_normal2);
    double signal_ = semi_normal1.dotProduct(r3);
    int signal;
    if (signal_ > 0){
        signal = 1;
    } else {
        signal = -1;
    };
    this->absValue = signal * angleD;
};

double Torsion::getValue(){
    return this->absValue;
};

vector <Point> Torsion::getPoints(){
    return vector <Point> {this->a, this->b, this->c, this->d};
};

void Torsion::increaseNdecrease(double increment, vector <char> freezePoints){
    Vector3D r = Vector3D(this->b.getCoords(), this->c.getCoords());
    vector <char> opt1 = {'a', 'b', 'c'};
    if (freezePoints == opt1){
        this->d.rotationVector(this->absValue+increment, r);
    } else {
        this->a.rotationVector(this->absValue+increment, r);
    };
    this->calcAbs();
};

