//
//  PeriodicTable.hh
//  MoleKing_util
//
//  Created by Thiago Lopes and Mateus Barbosa on 09/01/20.
//  Copyright © 2020 LMSC. All rights reserved.
//

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
