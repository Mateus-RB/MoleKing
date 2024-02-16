//   MoleKing //
//
//   File:        [Psi4Process.hpp]
//
//   Author(s):   ['LEEDMOL Reserch Group']
//   Site(s):     ['https://www.researchgate.net/lab/LEEDMOL-Heibbe-Cristhians-Research-Group-Heibbe-Cristhian-B-De-Oliveira']
//   Email(s):    ['heibbe@ufg.br']
//   Credits:     ['Copyright © 2023 LEEDMOL. All rights reserved.']
//   Date:        ['26.09.2023']
//   Version:     ['1.5.1']
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
    size_t beer;   
    size_t geometry; 
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
    string str_filePath;
    vector <string> geoStorage;
    vector <string> chargeStorage;
    vector <string> multiplicityStorage;


    //* molecule
    Molecule mol;

    //* bool

    bool ntFound;

    //* set functions
    void readOUTFile();
    void setMolecule();
    void setMul();
    void setCharge();

public:
    Psi4OUTfile(string filePath);   
    ~Psi4OUTfile();      
    string toStr();
    Molecule getMolecule();  
    int getMul();
    int getCharge();
};
