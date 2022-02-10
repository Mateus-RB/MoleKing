//
//  SupraMolecule.hh
//  MoleKing_util
//
//  Created by Thiago Lopes on 26/01/20.
//  Copyright © 2020 Laboratório de Estrutura Eletrônica e Dinâmica Molecular. All rights reserved.
//

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
    void setCharge();
    void getMoleculeBonds();
    void getMoleculeAngles();
    void getMoleculeTorsions();

public:
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
    void standardOrientation(int molNumber);
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
};


#endif /* SupraMolecule_hh */
