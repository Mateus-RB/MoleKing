//   CPP FILE HEADER //
//
//   File:        [G16Process.cpp]
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

#include "Testing.hpp"

/*
----------------------- G16LOGtest -----------------------
*/

G16LOGtest::G16LOGtest(string filePath, bool polarAsw){
    ifstream log_file(filePath);

    while (getline(log_file, line)){
        size_t scf = line.find("SCF Done:");
        size_t normalT = line.find("Normal termination of Gaussian");
        size_t basis = line.find("Standard basis:");
        size_t startMolecule = line.find("Input orientation:");
        size_t endMolecule = line.find("Distance matrix (angstroms):");
        
        if (scf != string::npos){
            value = line.substr(scf);

            size_t starterSCF = value.find("=");
            size_t endSCF = value.find("A.U.");

            size_t starterMethod = value.find("E(R");

            this->scf = stod(value.substr(starterSCF + 3, endSCF - starterSCF - 3));
            this->method = value.substr(starterMethod+3, starterSCF - starterMethod - 5);
        }
        
        if (basis != string::npos){
            this->basis = line.substr(basis + 16);
        };

        if (normalT != string::npos){
            this->info = line.substr(normalT + 31);
        }

        //make a string starting from startMolecule to endMolecule

        if (startMolecule != string::npos){
            string molecule = "";
            while (getline(log_file, line)){
                if (line.find("Distance matrix (angstroms):") != string::npos){
                    break;
                }
                molecule += line + "\n";
            }
            this->moleculeRange = molecule;
        }

    };
    log_file.close();
}

string G16LOGtest::getMol(){
    return this->moleculeRange;
};


/*
----------------------- G16LOGfile -----------------------
*/

string G16LOGtest::getDate(){
    return this->info;    
};

double G16LOGtest::getSCF(){
    return this->scf;
};

string G16LOGtest::getBasis(){
    return this->basis;
};

string G16LOGtest::getMethod(){
    return this->method;
};

string G16LOGtest::Summary(){
    return "Calculation done in " + this->info + "\nWith the level of theory " + this->method + "/" + this->basis + '.' + "\nSCF energy of " + to_string(this->scf) + " Hartree.";
};