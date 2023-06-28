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

    //* size_t
    size_t scf;
    size_t scfC;
    size_t homoFinder;
    size_t lumoFinder;
    size_t normalT;   
    size_t stdT; 
    size_t starterSCF;
    size_t endSCF;
    size_t starterMethod;  
    size_t basis; 

    //* string
    string basisValue;    
    string line;
    string value;
    string moleculeRange;
    string info;
    string method;
    string lumoSTR = "";
    string homoSTR = "";
    string moleculeSTR = "";
    string teste = "";
    string str_filePath;
    vector <string> iptStorage;
    vector <string> stdStorage;
    vector <string> homoStorage;
    vector <string> lumoStorage;

    //* double
    double scfValue;
    double homoValue;
    double lumoValue;

    //* molecule
    Molecule mol;
    

    //* Void
    void setMol();

    //* bool

    bool ntFound;
    bool stdFound;
    bool scfConvergence;

    //* map

    map<string, vector<string>> Orbitals;

    //* teste

   
public:
    G16LOGtest(string filePath, bool polarAsw = 0);     
    
    double getEnergy();
    double getHOMO(int index = -1);
    double getLUMO(int index = 0);
    map<string, vector<string>> getOrbitals(); 
    string getDate();    
    string getBasis();
    string getMethod();
    string getSummary();
    Molecule getMol();  
};