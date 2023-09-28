//   MoleKing //
//
//   File:        [G16LOGfile.hpp]
//
//   Author(s):   ['LEEDMOL Reserch Group']
//   Site(s):     ['https://www.researchgate.net/lab/LEEDMOL-Heibbe-Cristhians-Research-Group-Heibbe-Cristhian-B-De-Oliveira']
//   Email(s):    ['heibbe@ufg.br']
//   Credits:     ['Copyright © 2023 LEEDMOL. All rights reserved.']
//   Date:        ['17.01.2023']
//   Version:     ['1.5.0']
//   Status:      ['Development']
//   Language:    ['C++','Python']
//   Description: ['A python module written in C++ for theoretical chemestry']

#ifndef G16Process_hpp
#define G16Process_hpp

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

class G16LOGfile{
private:

    //* int

    int charge;
    int multiplicity;

    //* size_t
    size_t dipoleFinder;
    size_t mullikenFinder;
    size_t aHomoFinder;
    size_t aLumoFinder;
    size_t bHomoFinder;
    size_t bLumoFinder;
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
    vector <string> aHomoStorage;
    vector <string> aLumoStorage;
    vector <string> bHomoStorage;
    vector <string> bLumoStorage;
    vector <string> dipoleStorage;
    vector <string> tdStorage;

    vector<string> elstDipoleStorage;
    vector<string> alphaStorage;
    vector<string> betaStorage;
    vector<string> gammaStorage;

    //* double
    double scfValue;
    double aHomoValue;
    double aLumoValue;
    double bHomoValue;
    double bLumoValue;
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
    map<string, vector<string>> bOrbitals;
    map<int, map<string, double>> transitions;

    vector<string> Occupied;
    vector<string> bOccupied;
    vector<string> Unoccupied;
    vector<string> bUnoccupied;

    //* teste

    string aLumoStorageSTR = "";
    string aHomoStorageSTR = "";
    string bLumoStorageSTR = "";
    string bHomoStorageSTR = "";

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
    G16LOGfile(string filePath, bool polarAsw = 0, bool tdAsw = 0);   
    ~G16LOGfile();      
    map<string, vector<string>> getOrbitals(); 
    map<int, map<string, double>> getTransitions(const int index = 0);
    double getEnergy();
    vector<double> getHOMO(int index = -1);
    vector<double> getLUMO(int index = 0);    
    double getDipole(string axis = "tot");
    string getDate();    
    string getBasis();
    string getMethod();
    string toStr();
    Molecule getMolecule();  
};
