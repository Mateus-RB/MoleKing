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

//!----------------------- G16LOGtest -----------------------//

G16LOGtest::G16LOGtest(string filePath, bool polarAsw, bool tdAsw)
{
    // Initialize some variables
    ntFound = false;
    stdFound = false;
    scfConvergence = true;
    this->str_filePath = filePath;
    this->polarAsw = polarAsw;
    this->tdAsw = tdAsw;

    // SET FUNCTIONS
    readLOGFile();
    setMolecule();
    setDipole();
    setHOMO();
    setLUMO();
    setOrbitals();

    if (tdAsw)
    {
        setTransitions();
    }
};

//!----------------------- Set Functions -----------------------//

// Function to read file and extract the information

void G16LOGtest::readLOGFile()
{
    // Open the file at the given file path
    ifstream log_file(this->str_filePath);

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
        dipoleFinder = line.find("Tot=");
        tdFinder = line.find("Excited State ");

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

        // If the line contains "Convergence criterion not met.", set scfConvergence to false
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

        // Getting the dipoleValue
        if (dipoleFinder != string::npos)
        {
            // only append if "Dipole=" and "NTot=" not in line
            if ((line.find("Dipole=") == string::npos) && (line.find("NTot=") == string::npos))
            {
                this->dipoleStorage.push_back(line);
            }
        }

        // Getting HOMO and LUMO lines and appending into a vector;
        if (lumoFinder != string::npos)
        {
            this->lumoStorage.push_back(line);
        }

        if (homoFinder != string::npos)
        {
            this->homoStorage.push_back(line);
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
            this->info.pop_back(); // this remove the '.' of the string, wich is the last char in the string.
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

        // If the line contains "Excited State", extract the TD-DFT information
        if (tdAsw)
        {
            if (tdFinder != string::npos)
            {
                // only append if "Dipole=" and "NTot=" not in line
                if ((line.find("Excited State =") == string::npos) && (line.find("symmetry could not be determined.") == string::npos))
                {
                    this->tdStorage.push_back(line);
                }
            }
        }
    }

    // Close the file
    log_file.close();

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

    // If ntFound is false, throw an error
    if (!ntFound)
    {
        throw runtime_error("Normal termination of Gaussian not found in the log. Please check your log file.");
    }
}

// Function to set the molecule object using the extracted geometry
void G16LOGtest::setMolecule()
{
    stringstream ss(this->moleculeSTR);
    string line;

    while (getline(ss, line))
    {
        istringstream iss(line);
        vector<string> results((istream_iterator<string>(iss)), istream_iterator<string>());
        this->mol.addAtom(stoi(results[1]), stod(results[3]), stod(results[4]), stod(results[5]));
    }

    // If the molecule object has no atoms, throw an error
    if (this->mol.getSize() == 0)
    {
        throw runtime_error("Molecule not found in the log. Please check your output file.");
    }
};

void G16LOGtest::setOrbitals()
{
    this->Orbitals["Unoccupied"] = this->Unoccupied;
    this->Orbitals["Occupied"] = this->Occupied;
}

// Function to set the HOMO orbital of the calculation
void G16LOGtest::setHOMO()
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
    // this->Orbitals["Occupied"] = temp;
    this->Occupied = temp;
};

// Function to set the LUMO orbital of the calculation
void G16LOGtest::setLUMO()
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
    this->Unoccupied = temp;
};

// Function to set the transitions of the calculation
void G16LOGtest::setTransitions()
{
    if (!this->tdAsw)
    {
        throw runtime_error("Please, add tdAsw=1 to the G16LOGtest object to get the transitions.");
    }

    for (int i = 0; i < this->tdStorage.size(); i++)
    {
        stringstream ss(this->tdStorage[i]);
        string line;
        while (getline(ss, line))
        {
            if (line.find("Excited State") != string::npos)
            {
                istringstream iss(line);
                vector<string> results((istream_iterator<string>(iss)), istream_iterator<string>());

                string state = results[2];
                state.pop_back();
                int stateInt = stoi(state);
                string energy = results[4];
                double transitionEnergyDouble = stod(energy);
                string wv = results[6];
                double wvDouble = stod(wv);
                string osc = results[8];
                osc.erase(0, 2);
                double oscDouble = stod(osc);

                map<string, double> values = {
                    {"Energy", transitionEnergyDouble},
                    {"Wavelength", wvDouble},
                    {"Oscillation_Strength", oscDouble}};

                this->transitions[stateInt] = values;
            }
        }
    }
};

// Function to set the dipole moment of the calculation
void G16LOGtest::setDipole()
{
    stringstream ss(this->dipoleStorage.back());
    string line;

    while (getline(ss, line))
    {
        istringstream iss(line);
        vector<string> results((istream_iterator<string>(iss)), istream_iterator<string>());
        this->dipoleX = stod(results[1]);
        this->dipoleY = stod(results[3]);
        this->dipoleZ = stod(results[5]);
        this->dipoleTot = stod(results[7]);
    }
};

//!----------------------- Get Functions -----------------------//

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

// Function to user get the molecule object
Molecule G16LOGtest::getMolecule()
{
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

// Function to user get the orbitals
map<string, vector<string>> G16LOGtest::getOrbitals()
{
    return this->Orbitals;
};

// Function to user get the HOMO value;
double G16LOGtest::getHOMO(int index)
{
    if (index > 0)
    {
        // If the index is positive, throw an exception
        throw out_of_range("You've entered with a positive number. If you're looking for LUMO orbitals try to use: getLUMO(+n) function.");
    }

    else if (abs(index) == this->Occupied.size() || abs(index) > this->Occupied.size())
    {
        // If the index is greater than the number of HOMO orbitals found in the log file, throw an exception.
        throw out_of_range("Index out of range. You've entered with a number greater than the number of HOMO orbitals found in the log file. Number of orbital found: " + to_string(this->Occupied.size()) + ".");
    }

    else if (index == 0)
    {
        // If the index is 0, return the last value in the temp vector
        this->homoValue = stod(this->Occupied[this->Occupied.size() - 1]);
    }

    else
    {
        // subtract 1 from the index, the default one is zero then the index will be -1 [Highest occ molecular orbital], so if the user puts the -1, in C++ we'll look for the index -2 that will be the HOMO-1 orbital. Then return the value at that index in the temp vector
        index += -1;
        this->homoValue = stod(this->Occupied[this->Occupied.size() + index]);
    }
    return this->homoValue;
};

// Function to user get the LUMO value;
double G16LOGtest::getLUMO(int index)
{
    if (index < 0)
    {
        // If the index is negative, throw an exception
        throw out_of_range("You've entered with a negative number. If you're looking for HOMO orbitals try to use: getHOMO(-n) function.");
    }

    else if (index > this->Unoccupied.size() - 1)
    {
        // If the index is greater than the number of LUMO orbitals found in the log file, throw an exception
        throw out_of_range("Index out of range. You've entered with a number greater than the number of HOMO orbitals found in the log file. Number of orbital found: " + to_string(this->Unoccupied.size()) + ".");
    }

    else
    {
        // If the index is valid, return the value at that index in the temp vector
        this->lumoValue = stod(this->Unoccupied[index]);
    }

    return this->lumoValue;
};

// get Transitions
map<int, map<string, double>> G16LOGtest::getTransitions(int index)
{
    if (index < 0)
    {
        throw runtime_error("Do not exist Excited State lower than 1.");
    }
    else if (index > this->transitions.size())
    {
        throw runtime_error("You've entered with a number higher then the excited states found (" + to_string(this->transitions.size()) + ").");
    }

    else if (index == 0)
    {
        return this->transitions;
    }
    else
    {
        map<int, map<string, double>> temp;
        map<string, double> temp2;
        temp2 = {
            {"Energy", this->transitions[index]["Energy"]},
            {"Wavelength", this->transitions[index]["Wavelength"]},
            {"Oscillation_Strength", this->transitions[index]["Oscillation_Strength"]}};
        temp[index] = temp2;
        return temp;
    }
}

// Function to user get the Dipole Value
double G16LOGtest::getDipole(string axis)
{
    if (axis == "tot")
    {
        return this->dipoleTot;
    }

    else if (axis == "x")
    {
        return this->dipoleX;
    }

    else if (axis == "y")
    {
        return this->dipoleY;
    }

    else if (axis == "z")
    {
        return this->dipoleZ;
    }

    else
    {
        throw runtime_error("Invalid axis. Please, use 'tot', 'x', 'y' or 'z'.");
    }
};

string G16LOGtest::toStr()
{
    return "G16LOGFile: Calculation of " + this->str_filePath.substr(this->str_filePath.find_last_of("/") + 1) + " done in " + this->info + ", with the level of theory " + this->method + "/" + this->basisValue + " and " + "SCF energy of " + to_string(this->scfValue) + " Hartrees.";
}

//!----------------------- Destructor -----------------------//

G16LOGtest::~G16LOGtest()
{
    // clear the size_t
    this->dipoleFinder = 0;
    this->scf = 0;
    this->scfC = 0;
    this->homoFinder = 0;
    this->lumoFinder = 0;
    this->tdFinder = 0;
    this->normalT = 0;
    this->stdT = 0;
    this->starterSCF = 0;
    this->endSCF = 0;
    this->starterMethod = 0;
    this->basis = 0;

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
    this->charge = 0;

    // clear the molecules
    // this->mol.clear();

    // clear maps and vectors

    this->Orbitals.clear();
    this->Occupied.clear();
    this->Unoccupied.clear();
    this->transitions.clear();
};

//!----------------------- Notepad -----------------------//

//! TODO:

// TODO:getAlpha
// TODO:getBeta
// TODO:getGamma
// TODO:getOscillatorForce
// TODO:getWavelength
// TODO:getOscillatorForces
// TODO:getWavelengths
// TODO:getSymmetries
// TODO:getSymmetry
// TODO:getTransContributions
// TODO:getGradient
// TODO:getTransitions_str

//! Done:

//*  getEnergy | scfEnergy
//*  getMolecule
//*  getDipole
//*  getTransitions
//*  toStr()

//! New Features (didn't exist before):

//? getDate
//? getBasis
//? getMethod
//? getMolecule
//? getOrbitals
//? getHOMO() -> getHOMO(-n)
//? getLUMO -> getLUMO(+n)