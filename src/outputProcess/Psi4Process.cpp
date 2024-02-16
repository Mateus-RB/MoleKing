//   CPP FILE HEADER //
//
//   File:        [Psi4Process.cpp]
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


#include "Psi4Process.hpp"

//!----------------------- Psi4OUTfile -----------------------//

Psi4OUTfile::Psi4OUTfile(string filePath)
    : str_filePath(filePath)
{
    // SET FUNCTIONS
    readOUTFile();
    setCharge();
    setMul();
    setMolecule();
};

//!----------------------- Set Functions -----------------------//

// Function to read file and extract the information
void Psi4OUTfile::readOUTFile()
{
    // Open the file at the given file path
    ifstream log_file(this->str_filePath);
    // Loop through each line in the file
    while (getline(log_file, line))
    {
        // Check if the line contains certain keywords
        geometry = line.find("==> Geometry <==");
        chargeFinder = line.find("Charge       =");
        multiFinder = line.find("Multiplicity =");
        beer = line.find("Psi4 exiting successfully. Buy a developer a beer!");
        // Getting the charge and multiplicity
        if (chargeFinder != string::npos)
        {   
            this->chargeStorage.emplace_back(line);
        };
        if (multiFinder != string::npos)
        {   
            this->multiplicityStorage.emplace_back(line);
        };
        if (beer != string::npos)
        {
            ntFound = true;
        };

        if (line.find("==> Geometry <==") != string::npos)
        {
            while (getline(log_file, line))
            {
                if (line.find("Running in ") != string::npos)
                {
                    break;
                };
                moleculeSTR += line + "\n";
            };
            for (int i = 0; i < 8; i++)
            {
                moleculeSTR = moleculeSTR.substr(moleculeSTR.find("\n") + 1);
            };
            moleculeSTR = moleculeSTR.substr(0, moleculeSTR.rfind("\n"));
            geoStorage.emplace_back(moleculeSTR);
            moleculeSTR = "";
        };


    // Close the file
    };
    this->moleculeSTR = "";
    this->moleculeSTR = geoStorage[geoStorage.size() - 1];

    if (!ntFound)
    {
        throw runtime_error("Normal termination of Psi4 not found in the out. Please check your out file.");
    };

    log_file.close();
};

// Function to set the molecule object using the extracted geometry
void Psi4OUTfile::setMolecule()
{   
    stringstream ss(this->moleculeSTR);
    string line;
    int i = 0;
    while (getline(ss, line))
    {   
        istringstream iss(line);
        vector<string> results((istream_iterator<string>(iss)), istream_iterator<string>());
        this->mol.addAtom(results[0], stod(results[1]), stod(results[2]), stod(results[3]));
        i++;
    };
    this->mol.setCharge(this->charge);
    this->mol.setMultiplicity(this->multiplicity);
    // If the molecule object has no atoms, throw an error
    if (this->mol.getSize() == 0)
    {
        throw runtime_error("Molecule not found in the log. Please check your output file.");
    };
};

//!----------------------- Get Functions -----------------------//

// Function to user get the date and time of the calculation
void Psi4OUTfile::setMul()
{
    if (this->multiplicityStorage.size() > 0)
    {
        stringstream ss(this->multiplicityStorage.back());
        string line;
        while (getline(ss, line))
        {
            istringstream iss(line);
            vector<string> results((istream_iterator<string>(iss)), istream_iterator<string>());
            // extract the dipole moment components and convert them to the correct type
            this->multiplicity = stod(results[2]);

        };
    };
};

void Psi4OUTfile::setCharge()
{
    if (this->chargeStorage.size() > 0)
    {
        stringstream ss(this->chargeStorage.back());
        string line;
        while (getline(ss, line))
        {
            istringstream iss(line);
            vector<string> results((istream_iterator<string>(iss)), istream_iterator<string>());
            // extract the dipole moment components and convert them to the correct type
            this->charge = stod(results[2]);

        };
    };
};

int Psi4OUTfile::getMul()
{
    return this->multiplicity;
};

int Psi4OUTfile::getCharge()
{
    return this->charge;
};
Molecule Psi4OUTfile::getMolecule()
{
    return this->mol;
};


string Psi4OUTfile::toStr()
{
    return "Psi4OUTfile: Calculation of " + this->str_filePath.substr(this->str_filePath.find_last_of("/") + 1);
};

//!----------------------- Destructor -----------------------//

Psi4OUTfile::~Psi4OUTfile()
{
    //clear int
    this->charge = 0;
    this->multiplicity = 0;


    // clear the strings
    this->moleculeSTR.clear();

    // clear the molecules
    this->mol.clear();

};
