//
//  Hessian.hh
//  MoleKing_util
//
//  Created by Thiago Lopes and Mateus Barbosa on 21/01/20.
//  Copyright © 2020 LMSC. All rights reserved.
//

#ifndef Hessian_hpp
#define Hessian_hpp

#include <stdio.h>
#include "../math/Matrix.hpp"
#include "../math/Geometry.hpp"
#include <math.h>
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
