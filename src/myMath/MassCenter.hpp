//   MoleKing //
//
//   File:        [MassCenter.hpp]
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


#ifndef MassCenter_hpp
#define MassCenter_hpp
#include <stdio.h>
#include <vector>
#include "Geometry.hpp"

using namespace std;

class MassCenter{
    private:
    vector <double> massList;
    Point massCenterPoint;
    double axisMassCenter(vector <double> coords);

    public:
    MassCenter(vector <double> massList, vector <double> xCoords, vector <double> yCoords, vector <double> zCoords);
    ~MassCenter();
    Point getMassCenter();
};


#endif /* MassCenter_hpp */
