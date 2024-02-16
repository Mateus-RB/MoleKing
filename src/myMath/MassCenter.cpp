//   MoleKing //
//
//   File:        [MassCenter.cpp]
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

#include "MassCenter.hpp"

double MassCenter::axisMassCenter(vector <double> coords){
    double rUp = 0;
    double rDown = 0;
    for (int i = 0; i < (int) coords.size(); i++){
        rUp = rUp + (this->massList.at(i) * coords.at(i));
        rDown = rDown + this->massList.at(i);
    };
    return (rUp/rDown);
};

MassCenter::~MassCenter(){
    this->massList.clear();
    this->massCenterPoint.~Point();
    this->massList.resize(0);
};

MassCenter::MassCenter(vector <double> massList, vector <double> xCoords, vector <double> yCoords, vector <double> zCoords){
    this->massList = massList;
    this->massCenterPoint.setCoord('x', this->axisMassCenter(xCoords));
    this->massCenterPoint.setCoord('y', this->axisMassCenter(yCoords));
    this->massCenterPoint.setCoord('z', this->axisMassCenter(zCoords));
};

Point MassCenter::getMassCenter(){
        return this->massCenterPoint;
};

