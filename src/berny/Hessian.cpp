//   MoleKing //
//
//   File:        [Hessian.cpp]
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


#include "Hessian.hpp"

#define DAMPING 0.12
#define BOND_WEIGHT 0.45
#define ANGLE_WEIGHT 0.15
#define DIHEDRAL_WEIGHT 0.005

Hessian::Hessian(Molecule molecule){
    this->molecule = molecule;
}


double  Hessian::rho(vector <int> atoms){
    double bond = this->molecule.bondLength(atoms[0], atoms[1]);
    Atom atom1 = this->molecule.getAtomObj(atoms[0]);
    Atom atom2 = this->molecule.getAtomObj(atoms[1]);
    return  exp(- ( bond / (atom1.getAtomicRadio() + atom2.getAtomicRadio()) ) -1);
};

double  Hessian::wAngle(vector <int> atoms){
    double rho1 = this->rho(vector<int> {atoms[0], atoms[1]});
    double rho2 = this->rho(vector<int> {atoms[0], atoms[2]});
    return pow(rho1 * rho2, 0.5) * (DAMPING + (1 - DAMPING) * sin( this->molecule.valenceAngle(atoms[0], atoms[1], atoms[3]) ) );
};

double  Hessian::wDihedral(vector <int> atoms){
    double rho1 = this->rho(vector<int> {atoms[0], atoms[1]});
    double rho2 = this->rho(vector<int> {atoms[1], atoms[2]});
    double rho3 = this->rho(vector<int> {atoms[2], atoms[3]});
    return pow(rho1 * rho2 * rho3, 1/3) * (DAMPING + (1-DAMPING) * sin(this->molecule.valenceAngle(atoms[0], atoms[1], atoms[2]))) * (DAMPING + (1-DAMPING) * sin(this->molecule.valenceAngle(atoms[1], atoms[2], atoms[3])));
};


Matrix Hessian::doInitialGuess(){
    vector < vector <int> > bonds = molecule.getIRCBonds();
    vector < vector <int> > angles = molecule.getIRCAngles();
    vector < vector <int> > dihedrals = molecule.getIRCDihedrals();
    int dimension = (int) bonds.size() + (int) angles.size() + (int) dihedrals.size();
    Matrix hessian = Matrix(dimension, dimension);
    map < vector <int>, double > rhos;
    for (int i = 0; i < (int) bonds.size(); i++){
        rhos.insert(pair< vector <int>, double> ( bonds[i], this->rho(bonds[i]) ) );
    };
    vector <double> k_s(dimension);
    int a = 0;
    for (int i = 0; i < (int) bonds.size(); i++){
        k_s.at(a) = BOND_WEIGHT * rhos[bonds[i]];
        a += 1;
    };
    for (int i = 0; i < (int) angles.size(); i++){
        vector <int> ref1, ref2;
        if(angles[i][0] < angles[i][1]){
             ref1 = {angles[i][0], angles[i][1]};
        } else {
             ref1 = {angles[i][1], angles[i][0]};
        };
        if (angles[i][1] <  angles[i][2]){
            ref2 = {angles[i][1], angles[i][2]};
        } else {
            ref2 = {angles[i][2], angles[i][1]};
        };
        k_s.at(a) = ANGLE_WEIGHT * rhos[ref1] * rhos[ref2];
        a += 1;
    };
    for (int i = 0; i < (int) dihedrals.size(); i++){
        vector <int> ref1, ref2, ref3;
        if(dihedrals[i][0] < dihedrals[i][1]){
             ref1 = {dihedrals[i][0], dihedrals[i][1]};
        } else {
             ref1 = {dihedrals[i][1], dihedrals[i][0]};
        };
        if (dihedrals[i][1] <  dihedrals[i][2]){
            ref2 = {dihedrals[i][1], dihedrals[i][2]};
        } else {
            ref2 = {dihedrals[i][2], dihedrals[i][1]};
        };
        if(dihedrals[i][2] < dihedrals[i][3]){
             ref1 = {dihedrals[i][2], dihedrals[i][3]};
        } else {
             ref1 = {dihedrals[i][3], dihedrals[i][2]};
        };
        k_s.at(a) = DIHEDRAL_WEIGHT * rhos[ref1] * rhos[ref2];
        a += 1;
    };
    for (int i = 0; i < (int) k_s.size(); i++){
        hessian.replace(i, i, k_s[i]);
    };
    
    return hessian;
};

