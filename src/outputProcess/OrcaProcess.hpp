//   MoleKing //
//
//   File:        [G16LOGfile.hpp]
//
//   Author(s):   ['LEEDMOL Reserch Group']
//   Site(s):     ['https://www.researchgate.net/lab/LEEDMOL-Heibbe-Cristhians-Research-Group-Heibbe-Cristhian-B-De-Oliveira']
//   Email(s):    ['heibbe@ufg.br']
//   Credits:     ['Copyright © 2023 LEEDMOL. All rights reserved.']
//   Date:        ['17.01.2023']
//   Version:     ['1.6.0']
//   Status:      ['Development']
//   Language:    ['C++','Python']
//   Description: ['A python module written in C++ for theoretical chemistry']

#ifndef OrcaProcess_hpp
#define OrcaProcess_hpp

#include <stdio.h>
#include <string>
#include <fstream>
#include "../myMath/Vectors.hpp"
#include "../myMath/Matrix.hpp"
#include <iostream>
#include <numeric>
#include <regex>
#include "../chemicalUnits/Molecule.hpp"
#include <stdexcept>
#include <map>
#include <cmath>
#include "../chemicalUnits/PeriodicTable.hpp"

using namespace std;
#endif /* Testing_hpp */

class ORCALOGfile{
private:

    //* Constants

    PeriodicTable pt;

    //* istringstream

    istringstream logfile;

    //* int

    int charge;
    int multiplicity;

    //* size_t

    //* string

    string str_filePath;
    string line;
    string moleculeSTR = "";

    //* double

    double scfValue;
    double ZPE;
    double H;
    double G;
    double S;
    double sigma_r;
    double thetha_r;
    double qVib;
    double qRot;
    double qTrans;
    double qTot;
    double qEle;

    //* molecule

    Molecule mol;

    //* bool

    bool thermoAsw;
    bool isLinear;

    //* map and vectors

    vector <string> geometryBlocks;
    vector <string> vibFrequenciesStorage;
    vector <double> vibFrequencies;
    vector <double> vecThetha_r;

    //*test

    //* set functions

    void readLOGFile();
    void initializeLOGFile();
    void setMolecule();
    void setVibFrequencies();
    void setIsLinear();
    void set_qVib(double Temperature);
    void set_qRot(double Temperature);
    void set_qTrans(double Temperatura);
    void set_qTot();

public:
    ORCALOGfile(string filePath, bool thermoAsw = 0);
    ~ORCALOGfile();

    Molecule getMolecule();
    vector<double> getVibFrequencies();
    double getEnergy();
    double getZPE();
    double getH();
    double getG();
    double getS();
    double get_qEle();
    double get_qVib(double Temperature);
    double get_qRot(double Temperature);
    double get_qTrans(double Temperature);
    double get_qTot(double Temperature);
    double get_sigmaR();
};