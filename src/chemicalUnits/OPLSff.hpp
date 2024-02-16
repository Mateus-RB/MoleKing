//   MoleKing //
//
//   File:        [OPLSff.hpp]
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


#ifndef OPLSff_hpp
#define OPLSff_hpp

#include <stdio.h>
#include <map>
#include <iterator>
#include <string>
//#include "Molecule.hpp"

using namespace std;

class OPLSff{
    private:
    map<string, string> opls;
    map<string, double> epsilon;
    map<string, double> sigma;
    
    public:
    OPLSff();
    ~OPLSff();
    
    string getBondType(string opls_str);
    double getEpsilon(string opls_str);
    double getSigma(string opls_str);

};

#endif /* PeriodicTable_hpp */
