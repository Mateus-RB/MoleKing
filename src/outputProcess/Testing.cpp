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

//!----------------------- G16LOGtest -----------------------

G16LOGtest::G16LOGtest(string filePath, bool polarAsw)
{
    // Open the file at the given file path
    ifstream log_file(filePath);

    // Initialize some variables
    ntFound = false;
    stdFound = false;
    this->filePath = filePath;

    // Loop through each line in the file
    while (getline(log_file, line))
    {
        // Check if the line contains certain keywords
        scf = line.find("SCF Done:");
        normalT = line.find("Normal termination of Gaussian");
        stdT = line.find("Standard orientation:");
        basis = line.find("Standard basis:");

        // If the line contains "Normal termination of Gaussian", set ntFound to true
        if (normalT != string::npos)
        {
            ntFound = true;
        }

        // If the line contains "Standard orientation:", set stdFound to true
        if (stdT != string::npos)
        {
            stdFound = true;
        }

        // If the line contains "SCF Done:", extract the SCF value and calculation method used
        if (scf != string::npos)
        {
            value = line.substr(scf);

            starterSCF = value.find("=");
            endSCF = value.find("A.U.");

            starterMethod = value.find("E(R");

            this->scfValue = stod(value.substr(starterSCF + 3, endSCF - starterSCF - 3));
            this->method = value.substr(starterMethod + 3, starterSCF - starterMethod - 5);
        }

        // If the line contains "Standard basis:", extract the basis set used
        if (basis != string::npos)
        {
            this->basisValue = line.substr(basis + 16);
        };

        // If the line contains "Normal termination of Gaussian", extract the date and time of the calculation
        if (normalT != string::npos)
        {
            this->info = line.substr(normalT + 31);
        }

        // If the line contains "Input orientation:", extract the molecule's geometry in Input Orientation
        if (line.find("Input orientation:") != string::npos)
        {
            while (getline(log_file, line))
            {
                if (line.find("Distance matrix (angstroms):") != string::npos)
                {
                    break;
                }

                moleculeSTR += line + "\n";
            }

            for (int i = 0; i < 4; i++)
            {
                moleculeSTR = moleculeSTR.substr(moleculeSTR.find("\n") + 1);
            }

            moleculeSTR = moleculeSTR.substr(0, moleculeSTR.rfind("---------------------------------------------------------------------\n"));
            moleculeSTR = moleculeSTR.substr(0, moleculeSTR.rfind("\n"));

            iptStorage.push_back(moleculeSTR);

            moleculeSTR = "";
        }

        // If the line contains "Standard orientation:", extract the molecule's geometry in Standard Orientation
        if (line.find("Standard orientation:") != string::npos)
        {
            while (getline(log_file, line))
            {
                if (line.find(" Rotational constants (GHZ): ") != string::npos)
                {
                    break;
                }

                moleculeSTR += line + "\n";
            }

            for (int i = 0; i < 4; i++)
            {
                moleculeSTR = moleculeSTR.substr(moleculeSTR.find("\n") + 1);
            }

            moleculeSTR = moleculeSTR.substr(0, moleculeSTR.rfind("---------------------------------------------------------------------\n"));
            moleculeSTR = moleculeSTR.substr(0, moleculeSTR.rfind("\n"));

            stdStorage.push_back(moleculeSTR);

            moleculeSTR = "";
        }
    }

    // If stdFound is true, set the molecule string to the last geometry in stdStorage, as default molecule geometry is always be standard orientation if it is available
    if (stdFound)
    {   
        //TODO: When having more than one calculations like "code /home/phfmatias/Desktop/Testes/MoleKing/o_esch_Z_ts_Et_373.log", make moleking get all the tree geometries, not only the last one!

        this->moleculeSTR = "";
        this->moleculeSTR = stdStorage[stdStorage.size() - 1];
    }

    // If stdFound is false, set the molecule string to the last geometry in iptStorage
    else
    {
        this->moleculeSTR = "";
        this->moleculeSTR = iptStorage[iptStorage.size() - 1];
    }

    // Close the file
    log_file.close();

    // Set the molecule object using the extracted geometry
    setMol();

    // If ntFound is false, throw an error
    if (!ntFound)
    {
        throw runtime_error("Normal termination of Gaussian not found in the log. Please check your output file.");
    }
}

// Function to set the molecule object using the extracted geometry
void G16LOGtest::setMol()
{

    stringstream ss(this->moleculeSTR);
    string line;

    while (getline(ss, line))
    {
        istringstream iss(line);
        vector<string> results((istream_iterator<string>(iss)), istream_iterator<string>());
        this->mol.addAtom(stoi(results[1]), stod(results[3]), stod(results[4]), stod(results[5]));
    }
};

// Function to user get the date and time of the calculation
string G16LOGtest::getDate()
{
    return this->info;
};

// Function to user get the SCF energy
double G16LOGtest::getEnergy()
{
    return this->scfValue;
};

// Function to user get the basis set used
string G16LOGtest::getBasis()
{
    return this->basisValue;
};

// Function to user get the method used
string G16LOGtest::getMethod()
{
    return this->method;
};

// Function to user get a summary of the calculation
string G16LOGtest::getSummary()
{
    return "Calculation done in " + this->info + "\nWith the level of theory " + this->method + "/" + this->basisValue + '.' + "\nSCF energy of " + to_string(this->scfValue) + " Hartree.";
};

// Function to user get the molecule object
Molecule G16LOGtest::getMol()
{
    // If the molecule object has no atoms, throw an error
    if (this->mol.getSize() == 0)
    {
        throw runtime_error("Molecule not found in the log. Please check your output file.");
    }

    // If stdFound is true, print a message to the console
    if (this->stdFound)
    {   
        //check if this->filePath have / or not. if have, then split the string and get the last string, else just get the string
        string temp = this->filePath;
        if (temp.find("/") != string::npos)
        {
            temp = temp.substr(temp.find_last_of("/") + 1);
        }

        cout << "The geometry of " << temp << " was taken from the standard orientation." << endl;

    }

    return this->mol;
};