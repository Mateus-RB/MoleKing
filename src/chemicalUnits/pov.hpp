//   MoleKing //
//
//   File:        [pov.cpp]
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

#ifndef pov_hpp
#define pov_hpp

#include <stdio.h>
#include <string>
#include <fstream>
#include <iomanip>
#include <vector>
#include "Molecule.hpp"
#include "AtomicScale.hpp"
#include "PeriodicTable.hpp"

class PovRay{

private:
    void atomFarFromCM();
    string HEX2RGB(string hex);
    string getCamPos();
    string getLookAtPos();
    vector <double> setAtomicPosition(int atomicNumber);
    double setAtomicRadii(string atomicSymbol);
    string setAtomicColor(string atomicSymbol);
    Molecule mol;
    PeriodicTable PT;

    double maxDist;

    Point mc;
    vector <double> mcCoords;

public:

    PovRay(const Molecule& mol);
    void buildFile();
    
};

#endif /* pov_hh */