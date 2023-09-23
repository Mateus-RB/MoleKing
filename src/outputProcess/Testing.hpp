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

    //* int

    int charge;
    int multiplicity;

    //* size_t
    size_t dipoleFinder;
    size_t mullikenFinder;
    size_t homoFinder;
    size_t lumoFinder;
    size_t tdFinder;
    size_t scf;
    size_t scfC;
    size_t normalT;   
    size_t stdT; 
    size_t starterSCF;
    size_t endSCF;
    size_t starterMethod;  
    size_t basis; 
    size_t chargeMultiFinder;

    //* string
    string basisValue;    
    string line;
    string value;
    string moleculeRange;
    string info;
    string method;
    string moleculeSTR = "";
    string mullikenSTR = "";
    string str_filePath;
    vector <string> iptStorage;
    vector <string> stdStorage;
    vector <string> homoStorage;
    vector <string> lumoStorage;
    vector <string> dipoleStorage;
    vector <string> tdStorage;

    vector<string> elstDipoleStorage;
    vector<string> alphaStorage;
    vector<string> betaStorage;
    vector<string> gammaStorage;

    //* double
    double scfValue;
    double homoValue;
    double lumoValue;
    double dipoleTot;
    double dipoleX;
    double dipoleY;
    double dipoleZ;   

    //* molecule
    Molecule mol;

    //* bool

    bool ntFound;
    bool stdFound;
    bool scfConvergence;
    bool polarAsw;
    bool tdAsw;

    //* map and vectors

    map<string, vector<string>> Orbitals;
    map<int, map<string, double>> transitions;

    vector<string> Occupied;
    vector<string> Unoccupied;

    //* teste

    //* set functions
    void readLOGFile();
    void setMolecule();
    void setOrbitals();
    void setHOMO();
    void setLUMO();
    void setTransitions();
    void setDipole();
    void splitter();
    void setNLO();

public:
    G16LOGtest(string filePath, bool polarAsw = 0, bool tdAsw = 0);   
    ~G16LOGtest();      
    map<string, vector<string>> getOrbitals(); 
    map<int, map<string, double>> getTransitions(const int index = 0);
    double getEnergy();
    double getHOMO(int index = -1);
    double getLUMO(int index = 0);    
    double getDipole(string axis = "tot");
    string getDate();    
    string getBasis();
    string getMethod();
    string toStr();
    Molecule getMolecule();  
};
