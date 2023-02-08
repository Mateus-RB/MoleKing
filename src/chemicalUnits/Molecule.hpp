//   MoleKing //
//
//   File:        [Molecule.hpp]
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


#ifndef Molecule_hpp
#define Molecule_hpp

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
#include "../../include/Eigen/Eigenvalues"

class Molecule{

private:
    typedef vector<Atom> AtomList;
    double VDWRatio;
    AtomList molecule;
    typedef vector<ChargePoint> ChargeList;
    ChargeList chargePoint;
    typedef vector < vector <int> > VectorsInt;
    vector < pair < vector <int>, StraightSegment > > bonds;
    vector < pair < vector <int>, Angle > > angles;
    vector < pair < vector <int>, Torsion > > dihedrals;
    int multiplicity, charge;
    double angleToSpinInAref(int ref, char axisName);   
    void getBonds();
    void getAngles();
    void getDihedrals();
    vector<double> minNmaxValue(vector <double> v);
    vector<int> getBonded(int atomIndex);


public:
    typedef AtomList::iterator iterator;
    typedef AtomList::const_iterator const_iterator;

    Molecule();
    ~Molecule();


    Eigen::Matrix<double, 3, 3> moleculeTensor();
    void TesteTensor();

    void addChargePoints(double xPos, double yPos, double zPos, double charge);
    void addChargePoints(ChargePoint cp);
    void addAtom(string atomSymbol, double xPos, double yPos, double zPos, double atomicCharge = 0.0, bool freezeCode_ = 0);
    void addAtom(int atomNumber, double xPos, double yPos, double zPos, double atomicCharge = 0.0, bool freezeCode_ = 0);
    void addAtom(Atom atom);
    vector <string> getAtom(int number, bool symbol = 0);
    Atom getAtomObj(int number);
    ChargePoint getChargePointsObj(int number);
    void setCharge(int charge);
    double getCharge();
    long getSize();
    vector <Atom> getMoleculeVector();
    void setMultiplicity(int multiplicity);
    void setVDWRatio(double VDWratio);
    int getMultiplicity();
    double getVDWRatio();
    void normalizeCPs(int norm);
    vector< vector<string> > getMolecule(bool symbol = 0);
    vector< vector<string> > getChargePoints();
    Point getMassCenter();
    void spinMolecule(double angle, Vector3D spinVector);
    void spinMolecule(double angle, char axis);
    void translation(Vector3D traslationVector);
    void moveMassCenter(double x = 0.0, double y = 0.0, double z= 0.0);
    void moveTail(int atomNumber, double x = 0.0, double y = 0.0, double z= 0.0);
    void standardOrientation();
    vector <double> standardOrientationPath();
    double bondLength(int atomN1, int atomN2);
    double valenceAngle(int atomN1, int atomN2, int atomN3);
    double torsion(int atomN1, int atomN2, int atomN3, int atomN4);
    void doIRC();
    vector <int> molecularAxis();
    vector < vector <int> > getIRCBonds();
    vector < vector <int> > getIRCAngles();
    vector < vector <int> > getIRCDihedrals();
    Atom operator[](int index);
    iterator begin();
    iterator end();
    vector <Atom> moleculeList();
    void removeAtom(int atomNumber);
    void removeAtom(Atom atom);
    string toStr();
    bool operator==(Molecule mol);
    bool operator!=(Molecule mol);
    void removeElement(string element);
    Molecule copy();
    double getMolecularMass();
    //void doOPLS();
};

#endif /* Molecule_hh */
