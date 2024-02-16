//   MoleKing //
//
//   File:        [AtomicScale.hpp]
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


#ifndef AtomicScale_hpp
#define AtomicScale_hpp

#include "PeriodicTable.hpp"
#include "../myMath/Geometry.hpp"
#include <string>

class Atom{
    private:
    int atomicNumber;
    string atomicSymbol;
    double atomicMass;
    Point point;
    bool freezeCode;
    double atomicRadio;
    double charge;
    string opls;

    public:
    Atom(int atomicNumber, double x, double y, double z, double charge = 0.0, bool freezeCode_ = '0');
    Atom(string atomicSymbol, double x, double y, double z, double charge = 0.0, bool freezeCode_ = '0');
    double getAtomicMass();
    string getAtomicSymbol();
    int getAtomicNumber();
    bool operator==(Atom atom2);
    bool operator!=(Atom atom2);
    void setX(double newX);
    void setY(double newY);
    void setZ(double newZ);
    void setAtomicCharge(double newCharge);
    double getX();
    double getY();
    double getZ();
    double getAtomicCharge();
    double getAtomicRadio();
    void setNewPos(double newX, double newY, double newZ);
    void translation(Vector3D traslationVector);
    void rotationAxis(double tetha, Vector3D unitAxis);
    vector<double> getPos();
    bool operator<(Atom atom);
    bool operator>(Atom atom);
    bool operator<=(Atom atom);
    bool operator>=(Atom atom);
    Point getPoint();
    string toStr();
    int comp(Atom atom);
    void setOPLS(string opls_str);
    string getOPLS();
};

class ChargePoint{
    private:
    Point point;
    double charge;

    public:
    ChargePoint(double x, double y, double z, double charge);
    void setCharge(double newCharge);
    double getCharge();
    bool operator==(ChargePoint charge2);
    bool operator!=(ChargePoint charge2);
    void setX(double newX);
    void setY(double newY);
    void setZ(double newZ);
    double getX();
    double getY();
    double getZ();
    void setNewPos(double newX, double newY, double newZ);
    void translation(Vector3D traslationVector);
    void rotationAxis(double tetha, Vector3D unitAxis);
    vector<double> getPos();
    bool operator<(ChargePoint chargePoint);
    bool operator>(ChargePoint chargePoint);
    bool operator<=(ChargePoint chargePoint);
    bool operator>=(ChargePoint chargePoint);
    Point getPoint();
    string toStr();
    int comp(ChargePoint chargePoint);
};

#endif /* AtomicScale_hpp */
