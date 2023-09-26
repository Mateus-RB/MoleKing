//   MoleKing //
//
//   File:        [Psi4Process.hpp]
//
//   Author(s):   ['LEEDMOL Reserch Group']
//   Site(s):     ['https://www.researchgate.net/lab/LEEDMOL-Heibbe-Cristhians-Research-Group-Heibbe-Cristhian-B-De-Oliveira']
//   Email(s):    ['heibbe@ufg.br']
//   Credits:     ['Copyright © 2023 LEEDMOL. All rights reserved.']
//   Date:        ['26.09.2023']
//   Version:     ['1.5']
//   Status:      ['Development']
//   Language:    ['C++','Python']
//   Description: ['A python module written in C++ for theoretical chemestry']

#ifndef Psi4Process_hpp
#define Psi4Process_hpp

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

class Psi4OUTfile{
private:

    //* int

    int charge;
    int multiplicity;

    //* size_t
    size_t scf;
    size_t scfC;
    size_t beer;   
    size_t starterMethod;  
    size_t geo; 
    size_t multiFinder;
    size_t chargeFinder;

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
    vector <string> chargeStorage;
    vector <string> multiplicityStorage;
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
    void readOUTFile();
    void setMolecule();
    void splitter();
    void setMul();
    void setCharge();

public:
    Psi4OUTfile(string filePath);   
    ~Psi4OUTfile();      
    double getEnergy();
    string getDate();    
    string getBasis();
    string getMethod();
    string toStr();
    Molecule getMolecule();  
    int getMul();
    int getCharge();
};
