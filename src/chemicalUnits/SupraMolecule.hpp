//   MoleKing //
//
//   File:        [SupraMolecule.hpp]
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


#ifndef SupraMolecule_hpp
#define SupraMolecule_hpp

#include <stdio.h>
#include "Molecule.hpp"
#include <vector>


class SupraMolecule{

private:
    typedef vector<Molecule> MoleculeList;
    MoleculeList supraMolecule;
    typedef vector < vector <int> > VectorsInt;
    vector <VectorsInt> bonds;
    vector <VectorsInt> angles;
    vector <VectorsInt> dihedrals;
    int multiplicity, charge;
    void getMoleculeBonds();
    void getMoleculeAngles();
    void getMoleculeTorsions();
    void setCharge();

public:
    typedef MoleculeList::iterator iterator;
    typedef MoleculeList::const_iterator const_iterator;

    SupraMolecule();
    void addMolecule(Molecule molecule);
    void addAtomToMolecule(int molNumber, Atom atom);
    void addAtomToMolecule(int molNumber, string atomSymbol, double xPos, double yPos, double zPos);
    string toStr();
    Molecule getMolecule(int numberMolecule);
    void setMultiplicity(int multiplicity);
    int getMultiplicity();
    double getCharge();
    long getSize();
    Point getMassCenter();
    void translation(Vector3D traslationVector);
    void moveMassCenter(double x = 0.0, double y = 0.0, double z= 0.0);
    void moveTail(int molNumber, int atomNumber, double x = 0.0, double y = 0.0, double z = 0.0);
    void spinSupraMolecule(double angle, char axis);
    void spinSupraMolecule(double angle, Vector3D spinVector);
    //void standardOrientation(int molNumber);
    vector <VectorsInt> getIRCBonds();
    vector <VectorsInt> getIRCAngles();
    vector <VectorsInt> getIRCDihedrals();
    void removeAtom(int molNumber, int atomNumber);
    void removeAtom(Atom atom);
    void removeElement(int molNumber, string element);
    void removeElement(string element);
    void removeMolecule(int molNumber);
    void removeMolecule(Molecule molecule);
    double getSupraMolecularMass();
    vector <Molecule> moleculeList();
    iterator begin();
    iterator end();
    bool operator==(SupraMolecule sMol);
    bool operator!=(SupraMolecule sMol);
    bool operator<(SupraMolecule sMol);
    bool operator>(SupraMolecule sMol);
};


#endif /* SupraMolecule_hh */
