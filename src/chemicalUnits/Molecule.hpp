//   MoleKing //
//
//   File:        [Molecule.hpp]
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


#ifndef Molecule_hpp
#define Molecule_hpp

#include <stdio.h>
#include <string>
#include <fstream>
#include <iomanip>
#include "../myMath/MassCenter.hpp"
#include "AtomicScale.hpp"
#include "OPLSff.hpp"
#include "../myMath/Geometry.hpp"
#include "../myMath/Matrix.hpp"
#include <iterator>
#include <algorithm>
#include <vector>
#include "../../include/Eigen/Eigen"
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class Molecule{

private:
    Vector3D unitVector(Vector3D vector);
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
    void doIRC();
    string zmatrix = "";
    void reorderMolecule();


public:
    typedef AtomList::iterator iterator;
    typedef AtomList::const_iterator const_iterator;

    Molecule();
    ~Molecule();
    void toXYZ(string fileName = "MK_Molecule.xyz");
    void toGJF(string fileName = "MK_Molecule.gjf", string method = "B3LYP", string basis = "6-311g(d)", string addKeywords = "", string midKeywords = "", string endKeywords = "", int charge = 0, int multiplicity = 1, bool zmatrix=0, vector<double> EField = {});    
    void addChargePoints(double xPos, double yPos, double zPos, double charge);
    void addChargePoints(ChargePoint cp);
    void addAtom(string atomSymbol, double xPos, double yPos, double zPos, double atomicCharge = 0.0, bool freezeCode_ = 0);
    void addAtom(int atomNumber, double xPos, double yPos, double zPos, double atomicCharge = 0.0, bool freezeCode_ = 0);
    void addAtom(Atom atom);
    vector <string> getAtom(int number, bool symbol = 0);
    Atom& getAtomObj(int number);
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
    void moveMassCenter(double x = 0.0, double y = 0.0, double z = 0.0);
    void moveTail(int atomNumber, double x = 0.0, double y = 0.0, double z= 0.0);
    double bondLength(int atomN1, int atomN2);
    double valenceAngle(int atomN1, int atomN2, int atomN3);
    double torsion(int atomN1, int atomN2, int atomN3, int atomN4);
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
    bool operator<(Molecule mol);
    bool operator>(Molecule mol);
    void removeElement(string element);
    Molecule copy();
    double getMolecularMass();
    void clear();    
    double RMSD(Molecule MOL2);
    void alignMolecule(char axis = 'x');
    Eigen::Matrix3d getRotationMatrix(Eigen::Vector3d principalAxis, char mkAxis);

};

#endif /* Molecule_hh */
