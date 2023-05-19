//   MoleKing //
//
//   File:        [G16LOGtest.hpp]
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

#ifndef Testing_hpp
#define Testing_hpp

#include <stdio.h>
#include <string>
#include <fstream>
#include "../myMath/Vectors.hpp"
#include "../myMath/Matrix.hpp"
#include <iostream>
#include <regex>
#include "../chemicalUnits/Molecule.hpp"
#include <stdexcept>
#include <map>

using namespace std;
#endif /* Testing_hpp */

class G16LOGtest{
private:
    string line;
    string value;
    string moleculeRange;
    
    double scf;
    string info; 
    string basis;
    string method;
    

public:
    G16LOGtest(string filePath, bool polarAsw = 0);     
    double getSCF();
    string getDate();    
    string getBasis();
    string getMethod();
    string Summary();
    string getMol();   
};