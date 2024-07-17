//   MoleKing //
//
//   File:        [Hessian.hpp]
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


#ifndef Hessian_hpp
#define Hessian_hpp

#include <stdio.h>
#include "../myMath/Matrix.hpp"
#include "../myMath/Geometry.hpp"
#include <cmath>
#include <vector>
#include "../chemicalUnits/Molecule.hpp"
#include "../chemicalUnits/AtomicScale.hpp"
#include <iostream>

using namespace std;

class Hessian{
    
private:
    Matrix bondHessian, angleHessian, dihedralHessian;
    Molecule molecule;
    
public:
    Hessian(Molecule molecule);
    Matrix doInitialGuess();
    double rho(vector <int> atoms);
    double wAngle(vector <int> atoms);
    double wDihedral(vector <int> atoms);
    
};


#endif /* Hessian_hpp */
