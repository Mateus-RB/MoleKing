//   CPP FILE HEADER //
//
//   File:        [Psi4Process.cpp]
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


#include "Psi4Process.hpp"

//!----------------------- Psi4OUTfile -----------------------//

Psi4OUTfile::Psi4OUTfile(string filePath)
    : str_filePath(filePath)
{
    // SET FUNCTIONS
    readOUTFile();
    setCharge();
    setMul();
    //setMolecule();
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
        geo = line.find("==> Geometry <==");
        chargeFinder = line.find("Charge       =");
        multiFinder = line.find("Multiplicity =");      
        // Getting the charge and multiplicity
        if (chargeFinder != string::npos)
        {   
            this->chargeStorage.emplace_back(line);
        };
        if (multiFinder != string::npos)
        {   
            this->multiplicityStorage.emplace_back(line);
        };


    // Close the file
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
        this->mol.addAtom(stoi(results[1]), stod(results[3]), stod(results[4]), stod(results[5]));  
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
string Psi4OUTfile::getDate()
{
    return this->info;
};

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
// Function to user get the SCF energy 
double Psi4OUTfile::getEnergy()
{   
    // If scfConvergence is true, return the scfValue
    if (this->scfConvergence)
    {
        return this->scfValue;
    }
    // If scfConvergence is false, print a message to the console, to warn the user
    else
    {
        cerr << "WARNING in Psi4OUTfile::getEnergy(): SCF convergence not achieved. Please check your log file for 'Convergence criterion not met.'" << endl;
        return this->scfValue;
    };
};

// Function to user get the basis set used
string Psi4OUTfile::getBasis()
{
    return this->basisValue;
};

// Function to user get the method used
string Psi4OUTfile::getMethod()
{
    return this->method;
};

// Function to user get the molecule object
Molecule Psi4OUTfile::getMolecule()
{
    // If stdFound is true, print a message to the console
    if (!this->stdFound)
    {
        // check if this->filePath have / or not. if have, then split the string and get the last string, else just get the string
        string temp = this->str_filePath;
        if (temp.find("/") != string::npos)
        {
            temp = temp.substr(temp.find_last_of("/") + 1);
        };
        //cerr << "WARNING in Psi4OUTfile::getMolecule(): The geometry of " << temp << " was not taken from the standard orientation." << endl;
    };
    return this->mol;
};


string Psi4OUTfile::toStr()
{
    return "Psi4OUTfile: Calculation of " + this->str_filePath.substr(this->str_filePath.find_last_of("/") + 1) + " done in " + this->info + ", with the level of theory " + this->method + "/" + this->basisValue + " and " + "SCF energy of " + to_string(this->scfValue) + " Hartrees.";
};

//!----------------------- Destructor -----------------------//

Psi4OUTfile::~Psi4OUTfile()
{
    //clear int
    this->charge = 0;
    this->multiplicity = 0;

    // clear the size_t
    this->scf = 0;
    this->scfC = 0;
    this->starterMethod = 0;

    // clear the strings
    this->moleculeSTR.clear();

    // clear the doubless
    this->scfValue = 0;
    this->homoValue = 0;
    this->lumoValue = 0;
    this->dipoleTot = 0;
    this->dipoleX = 0;
    this->dipoleY = 0;
    this->dipoleZ = 0;

    // clear the molecules
    this->mol.clear();

    // clear maps and vectors

    this->Orbitals.clear();
    this->Occupied.clear();
    this->Unoccupied.clear();
    this->transitions.clear();
};
