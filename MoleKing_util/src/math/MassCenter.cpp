//
//  MassCenter.cpp
//  MoleKing_util
//
//  Created by Thiago Lopes and Mateus Barbosa on 09/01/20.
//  Copyright © 2020 LMSC. All rights reserved.
//

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

