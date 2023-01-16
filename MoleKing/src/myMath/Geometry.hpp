//
//  Geometry.hh
//  MoleKing_util
//
//  Created by Thiago Lopes and Mateus Barbosa on 10/01/20.
//  Copyright © 2020 LMSC. All rights reserved.
//

#ifndef Geometry_hpp
#define Geometry_hpp
#include <vector>
#include <math.h>
#include <iostream>
#include "Matrix.hpp"
#include "Vectors.hpp"

using namespace std;

class Point{
    private:
    double radius, tetha, phi, x, y, z;

    public:
    Point(double coord1, double coord2, double coord3, char typeCoord = 'c');
    Point();
    ~Point();
    bool operator==(Point point2);
    void setCoord(char coordName, double newValue);
    void setPoint(double coord1, double coord2, double coord3, char typeCoord = 'c');
    vector <double> getCoords(char typeCoord);
    void setCoords(vector <double> newValues, char typeCoord);
    void translation(Vector3D traslationVector);
    void rotationVector(double angle, Vector3D unitVector);
    string toStr(char spaceType = 'c');
};

class SphericalCoords{

    private:
    double x, y, z, radius, tetha, phi;

    public:
    SphericalCoords( double coord1/*x or radius*/, double coord2/*y or teta*/, double coord3/*z or phi*/, char spaceType);
    ~SphericalCoords();
    vector <double> toCartesian();
    vector <double> toSpherical();

};

class StraightSegment{
private:
    Point a, b;
    double absValue;
    void calcAbs();
    
public:
    StraightSegment(Point a, Point b);
    double getValue();
    vector <Point> getPoints();
    void stretchNcontract(double increment, char freezePoint = 'a');
};

class Angle{
private:
    Point a, b, c;
    double absValue;
    void calcAbs();
    
public:
    Angle(Point a, Point b, Point c);
    double getValue();
    vector <Point> getPoints();
    void increaseNdecrease(double increment, char freezePoint = 'a');
};

class Torsion{
private:
    Point a, b, c, d;
    double absValue;
    void calcAbs();
    
public:
    Torsion(Point a, Point b, Point c, Point d);
    double getValue();
    vector <Point> getPoints();
    void increaseNdecrease(double increment, vector <char> freezePoints = {'a', 'b'});
};

#endif /* Geometry_hpp */
