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
    scfConvergence = true;
    this->str_filePath = str_filePath;

    // Loop through each line in the file
    while (getline(log_file, line))
    {
        // Check if the line contains certain keywords
        scf = line.find("SCF Done:");
        normalT = line.find("Normal termination of Gaussian");
        stdT = line.find("Standard orientation:");
        scfC = line.find("Convergence criterion not met");
        basis = line.find("Standard basis:");
        homoFinder = line.find(" Alpha  occ. eigenvalues --");
        lumoFinder = line.find(" Alpha virt. eigenvalues --");

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

        // If the line contains ":", set scfConvergence to false

        if (scfC != string::npos)
        {
            this->scfConvergence = false;
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

        // Getting HOMO and LUMO lines and appending into a vector;

        if (lumoFinder != string::npos)
        {
            lumoSTR = "";
            lumoSTR += line;
            this->lumoStorage.push_back(lumoSTR);
        }

        if (homoFinder != string::npos)
        {
            homoSTR = "";
            homoSTR += line;
            this->homoStorage.push_back(homoSTR);
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
        this->moleculeSTR = "";
        this->moleculeSTR = stdStorage[stdStorage.size() - 1];
    }

    // If stdFound is false, set the molecule string to the last geometry in iptStorage
    else
    {
        this->moleculeSTR = "";
        this->moleculeSTR = iptStorage[iptStorage.size() - 1];
    }

    // Close the fileif (this->stdFound)
    {
        throw runtime_error("Normal termination of Gaussian not found in the log. Please check your log file.");
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

// Function to user get the HOMO value;

double G16LOGtest::getHOMO(int index)
{
    stringstream ss;
    vector<string> temp;

    for (int i = 0; i < this->homoStorage.size(); i++)
    {
        ss = stringstream(this->homoStorage[i]);
        string line;
        while (getline(ss, line))
        {
            istringstream iss(line);
            vector<string> results((istream_iterator<string>(iss)), istream_iterator<string>());
            for (int j = 0; j < results.size(); j++)
            {
                try
                {
                    // Check if the string can beif (this->stdFound)
                    // If the string cannot be converted to a double, ignore it
                }
            }
        }
    }

    if (index > 0)
    {
        // If the index is positive, throw an exception
        throw out_of_range("You entered with a positive number. If you're looking for LUMO orbitals try to use: getLUMO(+n) function.");
    }

    else if (abs(index) == temp.size() || abs(index) > temp.size())
    {
        // If the index is greater than the number of HOMO orbitals found in the log file, throw an exception.
        throw out_of_range("Index out of range. You entered with a number greater than the number of HOMO orbitals found in the log file. Number of orbital found: " + to_string(temp.size()));
    }

    else if (index == 0)
    {
        // If the index is 0, return the last value in the temp vector
        this->homoValue = stod(temp[temp.size() - 1]);
    }

    else
    {
        // subtract 1 from the index, the default one is zero then the index will be -1 [Highest occ molecular orbital], so if the user puts the -1, in C++ we'll look for the index -2 that will be the HOMO-1 orbital. Then return the value at that index in the temp vector
        index += -1;
        this->homoValue = stod(temp[temp.size() + index]);
    }

    return this->homoValue;
};

// Function to user get the LUMO value;

double G16LOGtest::getLUMO(int index)
{
    stringstream ss;
    vector<string> temp;

    for (int i = 0; i < this->lumoStorage.size(); i++)
    {
        ss = stringstream(this->lumoStorage[i]);
        string line;
        while (getline(ss, line))
        {
            istringstream iss(line);
            vector<string> results((istream_iterator<string>(iss)), istream_iterator<string>());
            for (int j = 0; j < results.size(); j++)
            {
                try
                {
                    // Check if the string can be converted to a double, then add to temp vector
                    stod(results[j]);
                    temp.push_back(results[j]);
                }
                catch (const std::exception &e)
                {
                    // If the string cannot be converted to a double, ignore it
                }
            }
        }
    }

    if (index < 0)
    {
        // If the index is negative, throw an exception
        throw out_of_range("You entered with a negative number. If you're looking for HOMO orbitals try to use: getHOMO(-n) function.");
    }

    else if (index > temp.size() - 1)
    {
        // If the index is greater than the number of LUMO orbitals found in the log file, throw an exception
        throw out_of_range("Index out of range. You entered with a number greater than the number of HOMO orbitals found in the log file. Number of orbital found: " + to_string(temp.size()));
    }

    else
    {
        // If the index is valid, return the value at that index in the temp vector
        this->lumoValue = stod(temp[index]);
    }

    return this->lumoValue;
};

// Function to user get the date and time of the calculation
string G16LOGtest::getDate()
{
    return this->info;
};

// Function to user get the SCF energy
double G16LOGtest::getEnergy()
{
    if (this->scfConvergence)
    {
        return this->scfValue;
    }

    else
    {
        cout << ("SCF convergence not achieved. Please, check your log file for -> 'Convergence criterion not met.'") << endl;
        return this->scfValue;
    }
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
    if (!this->stdFound)
    {
        // check if this->filePath have / or not. if have, then split the string and get the last string, else just get the string
        string temp = this->str_filePath;
        if (temp.find("/") != string::npos)
        {
            temp = temp.substr(temp.find_last_of("/") + 1);
        }

        cout << "The geometry of " << temp << " was not taken from the standard orientation." << endl;
    }

    return this->mol;
};