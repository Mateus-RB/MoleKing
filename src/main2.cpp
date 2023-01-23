//   MoleKing //
//
//   File:        [main.cpp]
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

#include <iostream>
#include <string>
#include <math.h>
#include <eigen3/Eigen/Eigenvalues>
#include "MoleKing.hpp"
#include "myMath/MassCenter.hpp"
#include "chemicalUnits/AtomicScale.hpp"
#include "myMath/Geometry.hpp"
#include "myMath/Matrix.hpp"
#include "chemicalUnits/Molecule.hpp"
#include "berny/Hessian.hpp"
#include "chemicalUnits/SupraMolecule.hpp"
#include "chemicalUnits/OPLSff.hpp"
#include "myMath/Vectors.hpp"
#include "outputProcess/G16Process.hpp"
using namespace std;

int main(){
    int AN1, AN2, AN3;
    double X1, X2, X3, Y1, Y2, Y3, Z1, Z2, Z3;
    cout<<"Add 3 Atoms" << endl;
    cout << "----------------" << endl;
    cout << "Atomic number of first Atom" << endl;
    cin >> AN1;
    cout << "X coordinated of first Atom" << endl;
    cin >> X1;
    cout << "Y coordinated of first Atom" << endl;
    cin >> Y1;
    cout << "Z coordinated of first Atom" << endl;
    cin >> Z1;
    cout << "----------------" << endl;
    cout << "Atomic number of second Atom" << endl;
    cin >> AN2;
    cout << "X coordinated of second Atom" << endl;
    cin >> X2;
    cout << "Y coordinated of second Atom" << endl;
    cin >> Y2;
    cout << "Z coordinated of second Atom" << endl;
    cin >> Z2;
    cout << "----------------" << endl;
    cout << "Atomic number of third Atom" << endl;
    cin >> AN3;
    cout << "X coordinated of third Atom" << endl;
    cin >> X3;
    cout << "Y coordinated of third Atom" << endl;
    cin >> Y3;
    cout << "Z coordinated of third Atom" << endl;
    cin >> Z3;
    Atom A1(AN1, X1, Y1, Z1);
    Atom A2(AN2, X2, Y2, Z2);
    Atom A3(AN3, X3, Y3, Z3);
    Molecule M = Molecule();
    M.addAtom(A1);
    M.addAtom(A2);
    M.addAtom(A3);
    cout << M.getVDWRatio() << endl;
    M.setVDWRatio(2.0);
    cout << M.getVDWRatio() << endl;
    cout << "The molecule is:" << endl;
    cout << M.toStr() << endl;
    
    return 0;  
};
