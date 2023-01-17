//
//  Molecule.hh
//  MoleKing_util
//
//  Created by Thiago Lopes and Mateus Barbosa on 21/01/20.
//  Copyright © 2020 LMSC. All rights reserved.
//
#ifndef Field_hpp
#define Field_hpp

#include <stdio.h>
#include <string>
#include "../myMath/MassCenter.hpp"
#include "AtomicScale.hpp"
#include "OPLSff.hpp"
#include "../myMath/Geometry.hpp"
#include "../myMath/Matrix.hpp"
#include <math.h>
#include <iterator>
#include <algorithm>
#include <vector>

class Field{

private:
    typedef vector<ChargePoint> ChargeList;
    ChargeList chargePoint;
    //typedef vector < vector <int> > VectorsInt;
    //vector < pair < vector <int>, StraightSegment > > bonds;
    //vector < pair < vector <int>, Angle > > angles;
    //vector < pair < vector <int>, Torsion > > dihedrals;
    //double angleToSpinInAref(int ref, char axisName);
    //void getBonds();
    //void getAngles();
    //void getDihedrals();
    //vector<double> minNmaxValue(vector <double> v);
    //vector<int> getBonded(int atomIndex);


public:
    typedef ChargeList::iterator iterator;
    typedef ChargeList::const_iterator const_iterator;

    Field();
    ~Field();
    void addChargePoints(double xPos, double yPos, double zPos, double charge);
    void addChargePoints(ChargePoint cp);
    //void addAtom(string atomSymbol, double xPos, double yPos, double zPos, double atomicCharge = 0.0, bool freezeCode_ = 0);
    //void addAtom(int atomNumber, double xPos, double yPos, double zPos, double atomicCharge = 0.0, bool freezeCode_ = 0);
    //void addAtom(Atom atom);
    //vector <string> getAtom(int number, bool symbol = 0);
    //Atom getAtomObj(int number);
    ChargePoint getChargePointsObj(int number);
    //void setCharge(int charge);
    //int getCharge();
    long getSize();
    //vector <Atom> getMoleculeVector();
    //void setMultiplicity(int multiplicity);
    //int getMultiplicity();
    void normalizeCPs(int norm);
    //vector< vector<string> > getMolecule(bool symbol = 0);
    vector< vector<string> > getChargePoints();
    //Point getMassCenter();
    void spinField(double angle, Vector3D spinVector);
    void spinFild(double angle, char axis);
    void translation(Vector3D traslationVector);
    void moveCenter(double x = 0.0, double y = 0.0, double z= 0.0);
    void moveTail(int indexCP, double x = 0.0, double y = 0.0, double z= 0.0);
    void standardOrientation();
    vector <double> standardOrientationPath();
    //double bondLength(int atomN1, int atomN2);
    //double valenceAngle(int atomN1, int atomN2, int atomN3);
    //double torsion(int atomN1, int atomN2, int atomN3, int atomN4);
    //void doIRC();
    //vector <int> molecularAxis();
    //vector < vector <int> > getIRCBonds();
    //vector < vector <int> > getIRCAngles();
    //vector < vector <int> > getIRCDihedrals();
    ChargePoint operator[](int index);
    iterator begin();
    iterator end();
    vector <ChargePoint> chargelist();
    void removeCP(int indexCP);
    void removeCP(ChargePoint cp);
    string toStr();
    bool operator==(Field field);
    bool operator!=(Field field);
    //void removeElement(string element);
    Field copy();
    //double getDipole();
    //void doOPLS();
};

#endif /* Molecule_hh */
