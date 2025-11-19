//   CPP FILE HEADER //
//
//   File:        [G16Process.cpp]
//
//   Author(s):   ['LEEDMOL Reserch Group']
//   Site(s):     ['https://www.researchgate.net/lab/LEEDMOL-Heibbe-Cristhians-Research-Group-Heibbe-Cristhian-B-De-Oliveira']
//   Email(s):    ['heibbe@ufg.br']
//   Credits:     ['Copyright © 2023 LEEDMOL. All rights reserved.']
//   Date:        ['17.01.2023']
//   Version:     ['1.5.4']
//   Status:      ['Development']
//   Language:    ['C++','Python']
//   Description: ['A python module written in C++ for theoretical chemestry']

#include "G16Process.hpp"

//!----------------------- G16LOGfile -----------------------//

G16LOGfile::G16LOGfile(string filePath, bool polarAsw, bool tdAsw, bool thermoAsw, bool cpAsw, double temperature, int link)
    : ntFound(false), stdFound(false), scfConvergence(true), str_filePath(filePath), polarAsw(polarAsw), tdAsw(tdAsw), thermoAsw(thermoAsw), cpAsw(cpAsw), temperature(temperature), link(link)
{
    // SET FUNCTIONS

    detectLink();

    if (this->linkStorage.size() == 1)
    {
       this->logfile.str(this->linkStorage[0]);
    }

    else if (this->linkStorage.size() > 1)
    {
        if (link > int(this->linkStorage.size()))
        {   
            throw runtime_error("ERROR in G16LOGfile::G16LOGfile(): Invalid log number. The number of link calculation(s) found in the log file is " + to_string(this->linkStorage.size()) + ".");
        }

        else if (link == 0)
        {
            //Throw an error if the user try to use link = 0. Tell him to use link = -1 to get the last link, or 1,2,3..
            throw runtime_error("ERROR in G16LOGfile::G16LOGfile(): Invalid log number. If you want to use the last link calculation, please use link = -1. To use the first link calculation, please use link = 1.");
        }

        else if (link == -1)
        {
            this->logfile.str(this->linkStorage[this->linkStorage.size() - 1]);
        }

        else if (link < 0 && link != -1)
        {
            throw runtime_error("ERROR in G16LOGfile::G16LOGfile(): Invalid log number. Log number can't be a negative number, except -1.");
        }

        else
        {
            this->logfile.str(this->linkStorage[link - 1]);
        };
    };

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

    if (polarAsw)
    {
        this->vecNLOFrec = {0.000};
        setNLO();
        setNLOFrequency();
        setAlpha();
        setBeta();
        setGamma();
    };

    if (thermoAsw)
    {
        setVibFrequencies();
        setIsLinear();
        setPrincipalAxesInertia();
        set_qVib();
        set_thetha_r();
        set_qRot();
        set_qTrans();
        set_qTot();
    };
};

//!----------------------- Set Functions -----------------------//


// Function to detect if the log have a link calculation

void G16LOGfile::detectLink()
{   

    normalT = line.find(" Normal termination of Gaussian");
    
    //check if the file exist

    ifstream log_file(this->str_filePath);

    if (!log_file.is_open())
    {
        throw runtime_error("ERROR in G16LOGfile::readLOGFile(): File not found. Please check your file path.");
    };

    while (getline(log_file, line))
    {
        linkStorageSTD += line + "\n";
        while (getline(log_file, line))
        {   
            if (line.find(" Normal termination") != string::npos)
            {
                linkStorageSTD += line + "\n";
                ntFound = true;
                break;
            };
            linkStorageSTD += line + "\n";
        };
        this->linkStorage.emplace_back(linkStorageSTD);
        linkStorageSTD = "";
    };    

    // If ntFound is false, throw an error
    if (!ntFound)
    {
        throw runtime_error("Normal termination of Gaussian not found in the log. Please check your log file.");
    };

    log_file.close();
}

// Function to read file and extract the information
void G16LOGfile::readLOGFile()
{
    while (getline(this->logfile, line))
    {   
        // Check if the line contains certain keywords
        scf = line.find("SCF Done:");
        stdT = line.find("Standard orientation:");
        scfC = line.find("Convergence criterion not met");
        basis = line.find("Standard basis:");
        aHomoFinder = line.find(" Alpha  occ. eigenvalues --");
        aLumoFinder = line.find(" Alpha virt. eigenvalues --");
        bHomoFinder = line.find(" Beta  occ. eigenvalues --");
        bLumoFinder = line.find(" Beta virt. eigenvalues --");
        dipoleFinder = line.find(" Tot=");
        cpFinder = line.find("Point Charges:");
        tdFinder = line.find("Excited State ");
        polarFinder = line.find(" Dipole moment:");
        chargeMultiFinder = line.find(" Charge =");
        normalT = line.find(" Normal termination of Gaussian");      

        // If the line contains "Standard orientation:", set stdFound to true
        if (stdT != string::npos)
        {
            stdFound = true;
        };
        // If the line contains "Convergence criterion not met.", set scfConvergence to false
        if (scfC != string::npos)
        {
            this->scfConvergence = false;
        };
        // If the line contains "SCF Done:", extract the SCF value and calculation method used
        if (scf != string::npos)
        {
            value = line.substr(scf);

            starterSCF = value.find("=");
            endSCF = value.find("A.U.");

            starterMethod = value.find("E(R");

            this->scfValue = stod(value.substr(starterSCF + 2, endSCF - starterSCF - 3));
            this->method = value.substr(starterMethod + 3, starterSCF - starterMethod - 5);
        };
        // Getting the charge and multiplicity
        if (chargeMultiFinder != string::npos)
        {   
            this->charge = stoi((line.substr(chargeMultiFinder + 11, 1)));
            this->multiplicity = stoi((line.substr(chargeMultiFinder + 28, 1)));
        };
        // Getting the dipoleValue
        if (dipoleFinder != string::npos)
        {
            // only append if "Dipole=" and "NTot=" not in line
            if ((line.find("Dipole=") == string::npos) && (line.find("NTot=") == string::npos))
            {
                this->dipoleStorage.emplace_back(line);
            };
        };
        // Getting HOMO and LUMO lines and appending into a vector;
        if (aHomoFinder != string::npos)
        {
            aHomoStorageSTR += line + "\n";
            while (getline(this->logfile, line))
            {   
                if (line.find(" Alpha virt. eigenvalues --") != string::npos)
                {
                    break;
                };
                aHomoStorageSTR += line + "\n";
            };
            this->aHomoStorage.emplace_back(aHomoStorageSTR);
            aHomoStorageSTR = "";
        };
        if (line.find(" Alpha virt. eigenvalues --") != string::npos) // Debugging :: aLumoFinder not match
        {
            aLumoStorageSTR += line + "\n";
            while (getline(this->logfile, line))
            {   
                if (line.find("          Condensed to atoms (all electrons):") != string::npos || line.find("     Molecular Orbital Coefficients:") != string::npos || line.find(" Electronic spatial extent (au):") != string::npos || line.find("  Beta  occ. eigenvalues --") != string::npos)
                {
                    bHomoStorageSTR += line + "\n";
                    break;
                };
                aLumoStorageSTR += line + "\n";
            };
            this->aLumoStorage.emplace_back(aLumoStorageSTR);
            aLumoStorageSTR = "";
        };

        if (bHomoFinder != string::npos)
        {
            bHomoStorageSTR += line + "\n";
            while (getline(this->logfile, line))
            {   
                if (line.find(" Beta virt. eigenvalues --") != string::npos)
                {
                    bLumoStorageSTR += line + "\n";
                    break;
                };
                bHomoStorageSTR += line + "\n";
            };
            this->bHomoStorage.emplace_back(bHomoStorageSTR);
            bHomoStorageSTR = "";
        };
        if (bLumoFinder != string::npos)
        {
            bLumoStorageSTR += line + "\n";
            while (getline(this->logfile, line))
            {   
                if (line.find("          Condensed to atoms (all electrons):") != string::npos || line.find("     Molecular Orbital Coefficients:") != string::npos || line.find(" Electronic spatial extent (au):") != string::npos)
                {
                    break;
                };
                bLumoStorageSTR += line + "\n";
            };
            this->bLumoStorage.emplace_back(bLumoStorageSTR);
            bLumoStorageSTR = "";
        };
        
        // If the line contains "Standard basis:", extract the basis set used
        if (basis != string::npos)
        {   
            this->basisValue = line.substr(basis + 16);
        };
        // If the line contains "Normal termination of Gaussian", extract the date and time of the calculation
        if (normalT != string::npos)
        {   
            this->info = line.substr(line.size() - 5);
            this->info.pop_back(); // this remove the '.' of the string, wich is the last char in the string.            
        };
        // If the line contains "Input orientation:", extract the molecule's geometry in Input Orientation
        //Teste!
        if (line.find("Input orientation:") != string::npos)
        {
            while (getline(this->logfile, line))
            {
                if (line.find("Distance matrix (angstroms):") != string::npos || line.find(" Rotational constants (GHZ): ") != string::npos || line.find("Stoichiometry") != string::npos || line.find(" Standard basis:") != string::npos)
                {
                    break;
                };
                moleculeSTR += line + "\n";
            };
            for (int i = 0; i < 4; i++)
            {
                moleculeSTR = moleculeSTR.substr(moleculeSTR.find("\n") + 1);
            };
            moleculeSTR = moleculeSTR.substr(0, moleculeSTR.rfind("---------------------------------------------------------------------\n"));
            moleculeSTR = moleculeSTR.substr(0, moleculeSTR.rfind(" ---------------------------------------------------------------------"));
            moleculeSTR = moleculeSTR.substr(0, moleculeSTR.rfind("\n"));
            iptStorage.emplace_back(moleculeSTR);
            moleculeSTR = "";
        };
        // If the line contains "Standard orientation:", extract the molecule's geometry in Standard Orientation
        if (line.find("Standard orientation:") != string::npos)
        {
            while (getline(this->logfile, line))
            {
                if (line.find("Distance matrix (angstroms):") != string::npos || line.find(" Rotational constants (GHZ): ") != string::npos || line.find(" Standard basis:") != string::npos)
                {
                    break;
                };
                moleculeSTR += line + "\n";
            };
            for (int i = 0; i < 4; i++)
            {
                moleculeSTR = moleculeSTR.substr(moleculeSTR.find("\n") + 1);
            };
            moleculeSTR = moleculeSTR.substr(0, moleculeSTR.rfind("---------------------------------------------------------------------\n"));
            moleculeSTR = moleculeSTR.substr(0, moleculeSTR.rfind(" ---------------------------------------------------------------------"));
            moleculeSTR = moleculeSTR.substr(0, moleculeSTR.rfind("\n"));
            stdStorage.emplace_back(moleculeSTR);
            moleculeSTR = "";
        };
        if (cpAsw)
        {
            if (cpFinder != string::npos)
            {
                //cpSTR += line + "\n";
                while (getline(this->logfile, line))
                {
                    if (line.find(" Pt Chg Charge= ") != string::npos)
                    {
                        break;
                    };
                    cpSTR += line + "\n";
                };
            };
        };

        // If the line contains "Excited State", extract the TD-DFT information
        if (tdAsw)
        {
            if (tdFinder != string::npos)
            {
                // only append if "Dipole=" and "NTot=" not in line
                if ((line.find("Excited State =") == string::npos) && (line.find("symmetry could not be determined.") == string::npos))
                {
                    this->tdStorage.emplace_back(line);
                };
            };
        };

        if (thermoAsw)
        {
            if (line.find(" Frequencies --") != string::npos)
            {
                vibFreqSTR += line + "\n";
                vibFrequenciesStorage.emplace_back(vibFreqSTR);
                vibFreqSTR = "";
            };
            
            if (line.find(" Principal axes and moments of inertia in atomic units:") != string::npos)
            {
                getline(this->logfile, line);
                getline(this->logfile, line);

                principalAxisSTR += line + "\n";
            };

            if (line.find("Sum of electronic and zero-point Energies=") != string::npos)
            {
                auto pos = line.find("=") + 1;
                string zpeString = line.substr(pos);
                zpeString.erase(0, zpeString.find_first_not_of(" \t"));
                zpeString.erase(zpeString.find_last_not_of(" \t\r\n") + 1);
                this->ZPE = stod(zpeString);
            };
            
            if (line.find(" Zero-point vibrational energy") != string::npos)
            {
                getline(this->logfile, line);
                // line ==                                    65.02685 (Kcal/Mol)
                auto pos = line.find(")") + 1;
                string zpveString = line.substr(0, pos);
                zpveString.erase(0, zpveString.find_first_not_of(" \t"));
                zpveString.erase(zpveString.find_last_not_of(" \t\r\n") + 1);
                this->ZPVE = stod(zpveString);
            };

            if (line.find("Rotational symmetry number") != string::npos)
            {
                // line ==  Rotational symmetry number  12.
                // get the last number before the '.'
                auto pos = line.find_last_of(" ");
                string sigmaString = line.substr(pos + 1);
                sigmaString.erase(sigmaString.find_last_not_of(".\t\r\n") + 1);
                this->sigma_r = stod(sigmaString);
            };

            if (line.find("Sum of electronic and thermal Enthalpies=") != string::npos)
            {
                auto pos = line.find("=") + 1;
                string enthalpyString = line.substr(pos);
                enthalpyString.erase(0, enthalpyString.find_first_not_of(" \t"));
                enthalpyString.erase(enthalpyString.find_last_not_of(" \t\r\n") + 1);
                this->H = stod(enthalpyString);
            };

            if (line.find("                     E (Thermal)             CV                S") != string::npos)
            {
                getline(this->logfile, line);
                getline(this->logfile, line);
                istringstream iss(line);
                vector<string> results((istream_iterator<string>(iss)), istream_iterator<string>());
                this->S = stod(results[results.size() - 1]);
            };

            if (line.find("Sum of electronic and thermal Free Energies=") != string::npos)
            {
                auto pos = line.find("=") + 1;
                string gibbsString = line.substr(pos);
                gibbsString.erase(0, gibbsString.find_first_not_of(" \t"));
                gibbsString.erase(gibbsString.find_last_not_of(" \t\r\n") + 1);
                this->G = stod(gibbsString);
            };

            if (line.find(" Electronic      0") != string::npos)
            {
                size_t dPos = line.find("D");
                if (dPos != string::npos)
                {
                    line.replace(dPos, 1, "e");
                };
                istringstream iss(line);
                vector<string> results((istream_iterator<string>(iss)), istream_iterator<string>());
                this->qEle = stod(results[1]);
            };
        };

        if (polarAsw)
        {   
            if (line.find("Electric dipole moment (input orientation):") != string::npos)
            {
                polarSTR += line + "\n";
                while (getline(this->logfile, line))
                {
                    if (line.find(" ----------------------------------------------------------------------") != string::npos)
                    {
                        break;
                    };
                    polarSTR += line + "\n";
                };
                polarStorageInp.emplace_back(polarSTR);
                polarSTR = "";
            };

            if (line.find(" Electric dipole moment (dipole orientation):") != string::npos)
            {
                this->dipFinder = true;
                polarSTR += line + "\n";
                while (getline(this->logfile, line))
                {
                    if (line.find(" ----------------------------------------------------------------------") != string::npos)
                    {
                        break;
                    };
                    polarSTR += line + "\n";
                };
                polarStorageDip.emplace_back(polarSTR);
                polarSTR = "";
            };
        };
    };

    // Close the file
    // log_file.close();

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
    };
};

// Function to set the molecule object using the extracted geometry
void G16LOGfile::setMolecule()
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
    // if cpAsw is true get the chargePoints from LOG
    if (this->cpAsw)
    {
        setChargePoints();
    };
};

void G16LOGfile::setChargePoints()
{   
    stringstream ss(this->cpSTR);
    string line;
    int i = 0;
    while (getline(ss, line))
    {   
        istringstream iss(line);
        vector<string> results((istream_iterator<string>(iss)), istream_iterator<string>());
        this->mol.addChargePoints(stod(results[1]), stod(results[2]), stod(results[3]), stod(results[5]));
        i++;
    };
};

void G16LOGfile::setOrbitals()
{   
    if (this->bOccupied.size() > 0 && this->bUnoccupied.size() > 0) {
        this->bOrbitals["Beta_Unoccupied"] = this->bUnoccupied;
        this->bOrbitals["Beta_Occupied"] = this->bOccupied;
        this->bOrbitals["Alpha_Unoccupied"] = this->Unoccupied;
        this->bOrbitals["Alpha_Occupied"] = this->Occupied;
    }
    else
    {
        this->Orbitals["Unoccupied"] = this->Unoccupied;
        this->Orbitals["Occupied"] = this->Occupied;
    }
};

// Function to set the HOMO orbital of the calculation
void G16LOGfile::setHOMO()
{
    vector<string> temp;
    string aHomoAuxiliary;
    aHomoAuxiliary = this->aHomoStorage[this->aHomoStorage.size() - 1];
    stringstream ss(aHomoAuxiliary);

    while (getline(ss, line))
    {   
        istringstream iss(line);
        vector<string> results ((istream_iterator<string>(iss)), istream_iterator<string>());

        for (int i = 4; i < results.size(); i++)
        {   
            try
            {   
                temp.emplace_back(results[i]);
            }
            catch (const std::exception &e)
            {
                // If the string cannot be converted to a double, ignore it
            };
        }
    };
    
    if (this->bHomoStorage.size() > 0)
    {
        vector<string> bTemp;
        string bHomoAuxiliary;
        bHomoAuxiliary = this->bHomoStorage[this->bHomoStorage.size()-1];
        stringstream ss(bHomoAuxiliary);

        while (getline(ss, line))
        {   
            istringstream iss(line);
            vector<string> bResults ((istream_iterator<string>(iss)), istream_iterator<string>());

            for (int i = 4; i < bResults.size(); i++)
            {   
                try
                {   
                    bTemp.emplace_back(bResults[i]);
                }
                catch (const std::exception &e)
                {
                    // If the string cannot be converted to a double, ignore it
                };
            }
        };  
        this->bOccupied = bTemp;
    }
    this->Occupied = temp;
};

// Function to set the LUMO orbital of the calculation
void G16LOGfile::setLUMO()
{   
    vector<string> temp;
    string aLumoAuxiliary;
    aLumoAuxiliary = this->aLumoStorage[this->aLumoStorage.size() - 1];
    stringstream ss(aLumoAuxiliary);

    while (getline(ss, line))
    {   
        istringstream iss(line);
        vector<string> results ((istream_iterator<string>(iss)), istream_iterator<string>());

        for (int i = 4; i < results.size(); i++)
        {   
            try
            {   

                temp.emplace_back(results[i]);
            }
            catch (const std::exception &e)
            {
                // If the string cannot be converted to a double, ignore it
            };
        }
    }

    if (this->bLumoStorage.size() > 0)
    {
        vector<string> bTemp;
        string bLumoAuxiliary;
        bLumoAuxiliary = this->bLumoStorage[this->bLumoStorage.size()-1];
        stringstream ss(bLumoAuxiliary);

        while (getline(ss, line))
        {   
            istringstream iss(line);
            vector<string> bResults ((istream_iterator<string>(iss)), istream_iterator<string>());

            for (int i = 4; i < bResults.size(); i++)
            {   
                try
                {   
                    bTemp.emplace_back(bResults[i]);
                }
                catch (const std::exception &e)
                {
                    // If the string cannot be converted to a double, ignore it
                };
            }
        };  
        this->bUnoccupied = bTemp;
    }

    this->Unoccupied = temp;
};

// Function to set the transitions of the calculation
void G16LOGfile::setTransitions()
{   
    //if tdAsw is false, throw an error
    if (!this->tdAsw)
    {
        throw runtime_error("Please, add tdAsw=1 to the G16LOGfile object to get the transitions.");
    }
    // if true, iterate through the tdStorage vector
    for (int i = 0; i < this->tdStorage.size(); i++)
    {
        stringstream ss(this->tdStorage[i]);
        string line;
        while (getline(ss, line))
        {   
            // if the line contains "Excited State", extract the information
            if (line.find("Excited State") != string::npos)
            {
                istringstream iss(line);
                vector<string> results((istream_iterator<string>(iss)), istream_iterator<string>());
                // extract the state, energy, wavelength and oscillator strength, and conver it to the correct type
                string state = results[2];
                state.pop_back();
                int stateInt = stoi(state);
                string energy = results[4];
                double transitionEnergyDouble = stod(energy);
                string wv = results[6];
                double wvDouble = stod(wv);
                string osc = results[8];

                //this is removing the first two characters of the string, which are the F=
                osc.erase(0, 2);
                double oscDouble = stod(osc);

                // create a map with the extracted information and add it to the transitions map
                map<string, double> values = {
                    {"Energy", transitionEnergyDouble},
                    {"Wavelength", wvDouble},
                    {"Oscillation_Strength", oscDouble}};
                this->transitions[stateInt] = values;
            };
        };
    };
};

// Function to set the dipole moment of the calculation
void G16LOGfile::setDipole()
{   
    if (this->dipoleStorage.size() > 0)
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
        };
    };
};

//!----------------------- Get Functions -----------------------//

// Function to user get the date and time of the calculation
string G16LOGfile::getDate()
{
    return this->info;
};

// Function to user get the SCF energy 
double G16LOGfile::getEnergy()
{   
    // If scfConvergence is true, return the scfValue
    if (this->scfConvergence)
    {
        return this->scfValue;
    }
    // If scfConvergence is false, print a message to the console, to warn the user
    else
    {
        cerr << "WARNING in G16LOGfile::getEnergy(" << this->str_filePath << "): SCF convergence not achieved. Please check your log file for 'Convergence criterion not met.'" << endl;
        return this->scfValue;
    };
};

// Function to user get the basis set used
string G16LOGfile::getBasis()
{
    return this->basisValue;
};

// Function to user get the method used
string G16LOGfile::getMethod()
{
    return this->method;
};

// Function to user get the molecule object
Molecule G16LOGfile::getMolecule()
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
    };
    return this->mol;
};

// Function to user get the orbitals
map<string, vector<string>> G16LOGfile::getOrbitals() 
{   
    if (this->bOccupied.size() > 0 && this->bUnoccupied.size() > 0)
    {   
        //throw an warning
        return this->bOrbitals;
    }
    else
    {
        return this->Orbitals;

    }
};

// Function to user get the HOMO value;
vector<double> G16LOGfile::getHOMO(int index) 
{   
    if (index > 0)
    {
        // If the index is positive, throw an exception
        throw out_of_range("ERROR in G16LOGfile::getHOMO(): Invalid index. Index can't be a positive number. If you're looking for LUMO orbitals try to use: getLUMO(+n) function.");
    }

    else if (abs(index) == this->Occupied.size() || abs(index) > this->Occupied.size())
    {
        // If the index is greater than the number of HOMO orbitals found in the log file, throw an exception.
        throw out_of_range("ERROR in G16LOGfile::getHOMO(): Index out of range. Index greather than the number greater than the number of HOMO orbitals found in the log file. Number of orbital found: " + to_string(this->Occupied.size()) + ".");
    }

    else if (index == 0)
    {
        // If the index is 0, return the last value in the temp vector
        if (this->bOccupied.size() > 0)
        {
            this->bHomoValue = stod(this->bOccupied[this->bOccupied.size() - 1]);
        }
        this->aHomoValue = stod(this->Occupied[this->Occupied.size() - 1]);
    }

    else
    {
        // subtract 1 from the index, the default one is zero then the index will be -1 [Highest occ molecular orbital], so if the user puts the -1, in C++ we'll look for the index -2 that will be the HOMO-1 orbital. Then return the value at that index in the temp vector
        index += -1;

        if (this->bOccupied.size() > 0)
        {
            this->bHomoValue = stod(this->bOccupied[this->bOccupied.size() + index]);
        }
        
        this->aHomoValue = stod(this->Occupied[this->Occupied.size() + index]);
    };

    if (this->bOccupied.size() > 0)
    {   
        cerr << "WARNING in G16LOGfile::getHOMO(). Your job is an Open Shell. Please, use as alpha_HOMO, beta_HOMO = getHOMO()." << endl;
        return {this->aHomoValue, this->bHomoValue};
    } 

    else 
    {
        return {this->aHomoValue};
    }
};

// Function to user get the LUMO value;
vector<double> G16LOGfile::getLUMO(int index)
{
    if (index < 0)
    {
        // If the index is negative, throw an exception
        throw out_of_range("ERROR in G16LOGfile::getLUMO(): Invalid index. Index can't be a negative number. If you're looking for HOMO orbitals try to use: getHOMO(-n) function.");
    }

    else if (index > this->Unoccupied.size() - 1)
    {
        // If the index is greater than the number of LUMO orbitals found in the log file, throw an exception
        throw out_of_range("ERROR in G16LOGfile::getLUMO(): Index out of range.  Index greater than the number of HOMO orbitals found in the log file. Number of orbital found: " + to_string(this->Unoccupied.size()) + ".");
    }

    else
    {
        // If the index is valid, return the value at that index in the temp vector
        if (this->bUnoccupied.size() > 0)
        {
            this->bLumoValue = stod(this->bUnoccupied[index]);
        }
        this->aLumoValue = stod(this->Unoccupied[index]);
    };

    if (this->bUnoccupied.size() > 0)
    {   
        cerr << "WARNING in G16LOGfile::getLUMO(). Your job is an Open Shell. Please, use as alpha_LUMO, beta_LUMO = getLUMO()." << endl;
        return {this->aLumoValue, this->bLumoValue};
    } 

    else 
    {
        return {this->aLumoValue};
    }
};

// get Transitions
map<int, map<string, double>> G16LOGfile::getTransitions(int index)
{
    if (index < 0)
    {
        throw runtime_error("ERROR in G16LOGfile::getTransitions(): Invalid index. Excited state indices start at 1.");
    }
    else if (index > this->transitions.size())
    {
        throw runtime_error("ERROR in G16LOGfile::getTransitions(): Invalid index. The number of excited states found in the log file is " + to_string(this->transitions.size()) + ".");
    }

    else if (index == 0)
    {
        if (this->tdAsw && this->transitions.size() == 0)
        {
            throw runtime_error("ERROR in G16LOGfile::getTransitions(): No transitions found in the log file.");
        };

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

        if (this->tdAsw && temp.size() == 0)
        {
            throw runtime_error("ERROR in G16LOGfile::getTransitions(): No transitions found in the log file.");
        };

        return temp;
    };
};

double G16LOGfile::getDipole(string axis)
{   
    if (this->dipoleStorage.size() == 0)
    {
        throw runtime_error("ERROR in G16LOGfile::getDipole(): No dipole found in the log file.");
    }    

    else if (axis == "tot")
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
        throw runtime_error("ERROR in G16LOGfile::getDipole(): Invalid axis. Please, use 'tot', 'x', 'y' or 'z'.");
    };
};

// Function to user get the frequency of the calculation

void G16LOGfile::setNLO()
{   
    if ( this->polarStorageDip.size() == 0 )
    {
        cerr << "WARNING in G16LOGfile: Your molecule is symmetric. So your dipole is zero, and there is no dipole orientation." << endl;
    };

    if ( (this->polarStorageDip.size() == 0) && (this->polarStorageInp.size() == 0) )
    {
        throw runtime_error("ERROR in G16LOGfile::setNLO(): No NLO found in the log file.");
    };
    
    if (this->polarStorageDip.size() > 0) 
    {
        // get the last element of the polarStorage vector
        this->polarAuxiliaryDip = this->polarStorageDip[this->polarStorageDip.size() - 1];
        // create a stringstream object with polarAuxiliary
        stringstream sdip(polarAuxiliaryDip);
        // Store each line of the stringstream in the lines vector
        while (getline(sdip, line))
        {
            this->vecPolarDip.emplace_back(line);
        };
        //check if the last element of the vector lines is equal to "" or " " or "\n"
        if (this->vecPolarDip[vecPolarDip.size() - 1] == "" || this->vecPolarDip[vecPolarDip.size() - 1] == " " || this->vecPolarDip[vecPolarDip.size() - 1] == "\n")
        {
            // if true, remove the last element of the vector lines
            this->vecPolarDip.pop_back();
        };
    };
    if (this->polarStorageInp.size() == 0){
        throw runtime_error("ERROR in G16LOGfile::setNLO(): No NLO found in the log file.");
    }
    if (this->polarStorageInp.size() > 0)
    {   
        this->polarAuxiliaryInp = this->polarStorageInp[this->polarStorageInp.size() - 1];
        stringstream sinp(polarAuxiliaryInp);
        
        while (getline(sinp, line))
        {
            this->vecPolarInp.emplace_back(line);
        };
        if (this->vecPolarInp[vecPolarInp.size() - 1] == "" || this->vecPolarInp[vecPolarInp.size() - 1] == " " || this->vecPolarInp[vecPolarInp.size() - 1] == "\n")
        {
            this->vecPolarInp.pop_back();

        };
    }
};

vector<string> G16LOGfile::getNLO(string Orientation)
{   
    if (!this->polarAsw)
    {
        throw runtime_error("ERROR in G16LOGfile::getNLO(): No NLO found in the log file., Try using polarAsw=1.");
    };
    if(Orientation == "input")
    {
        return this->polarStorageInp;
    }
    
    if(Orientation == "dipole")
    {   
        return this->polarStorageDip;
    }
    else
    {
        throw runtime_error("ERROR in G16LOGfile::getNLO(): Invalid Orientation. Please, use 'input' or 'dipole'.");
    };
};

void G16LOGfile::setNLOFrequency()
{    
    if (this->dipFinder)
    {
        for (int i = 0; i < this->vecPolarDip.size(); i++)
        {
            if (this->vecPolarDip[i].find("Alpha(-w;w)") != string::npos)
            {
                this->Freq = this->vecPolarDip[i].substr(this->vecPolarDip[i].find("Alpha(-w;w)") + 14);
                this->Freq = this->Freq.substr(0, this->Freq.find("nm"));
                //convert Freq to double
                this->FreqDouble = stod(this->Freq);
                // add the frequency to the vecNLOFrec vector
                this->vecNLOFrec.emplace_back(this->FreqDouble);
            };
        };
    }

    else
    {
        for (int i = 0; i < this->vecPolarInp.size(); i++)
        {
            if (this->vecPolarInp[i].find("Alpha(-w;w)") != string::npos)
            {
                this->Freq = this->vecPolarInp[i].substr(this->vecPolarInp[i].find("Alpha(-w;w)") + 14);
                this->Freq = this->Freq.substr(0, this->Freq.find("nm"));
                //convert Freq to double
                this->FreqDouble = stod(this->Freq);
                // add the frequency to the vecNLOFrec vector
                this->vecNLOFrec.emplace_back(this->FreqDouble);
            };
        };
    };
};

vector<double> G16LOGfile::getNLOFrequency()
{   
    if (!this->polarAsw)
    {
        throw runtime_error("ERROR in G16LOGfile::getNLOFrequency(): No frequency found in the log file. Try using polarAsw=1.");
    }
    return this->vecNLOFrec;
    
};

void G16LOGfile::setAlpha()
{   
    string WhatPolar = "";
    vector<string> UsePolar;

    for( int o = 0; o < 2; o++){
        map<double, int> start, end;
        map<double,map<string,vector<string>>> FrequencyInfo;
        map<string,vector<string>> NLOInfo;
        vector<vector<string>> Alpha_0, Alpha_w;

        if (o==0){
            WhatPolar = "input";
            UsePolar = this->vecPolarInp;
        }
        else
        {
            WhatPolar = "dipole";
            UsePolar = this->vecPolarDip;
        }

        for (int f = 0; f < this->vecNLOFrec.size(); f++){
            double freqPrincipal = this->vecNLOFrec[f];

            for (int i = 0; i < UsePolar.size(); i++)
            {
                if (UsePolar[i].find("Alpha(0;0)") != string::npos)
                {
                    start.insert(make_pair(this->vecNLOFrec[0],i+2));
                    end.insert(make_pair(this->vecNLOFrec[0],i+10));
                }
                else if (UsePolar[i].find("Alpha(-w;w) w= ") != string::npos)
                {   
                    // OK
                    string Freq = UsePolar[i].substr(UsePolar[i].find("Alpha(-w;w)") + 14);
                    Freq = Freq.substr(0, Freq.find("nm"));
                    //convert Freq to double
                    double FreqDouble = stod(Freq);
                    if(freqPrincipal == FreqDouble)
                    {
                        start.insert(make_pair(FreqDouble,i+2));
                        end.insert(make_pair(FreqDouble,i+10));
                    }
                }
            }
        for (auto it = start.begin(); it != start.end(); ++it) 
        {   
            for(int i = start[it->first]; i < end[it->first]; i++)
            {   
                vector<string> T = this->customSplit(UsePolar[i]);         
                NLOInfo.insert(make_pair(T[0], vector<string>({T[1], T[2], T[3]})));
            };
            FrequencyInfo.insert(make_pair(it->first, NLOInfo));
            NLOInfo.clear();
        }    
    };
    this->Alpha.insert(make_pair(WhatPolar,FrequencyInfo));
    };
};


map<string,double> G16LOGfile::getAlpha(string orientation, string unit, double frequency)
{   
    map<string,double> temp;
    map<double,map<string,vector<string>>> alpha;

    std::transform(orientation.begin(), orientation.end(), orientation.begin(),[](unsigned char c){ return std::tolower(c); });
        
    if (orientation == "input"){
        alpha = this->Alpha[orientation];
    }
    else{
        alpha = this->Alpha["dipole"];
    }
    if (!this->polarAsw)
    {
        throw runtime_error("ERROR in G16LOGfile::getAlpha(): No NLO found in the log file., Try using polarAsw=1.");
    }

    //check if frequency is in the vecNLOFrec vector
    if (find(this->vecNLOFrec.begin(), this->vecNLOFrec.end(), frequency) == this->vecNLOFrec.end())
    {   
        string temp = "";

        for (int i = 0; i < this->vecNLOFrec.size(); i++)
        {
            temp += to_string(this->vecNLOFrec[i]) + ", ";
        }
        
        throw runtime_error("ERROR in G16LOGfile::getAlpha(): Frequency not found in the log file. Try: " + temp + "instead.");
    }
    
    //get the element with the frequency key from map this->Alpha
    //this->Alpha[frequency];

    //create a new map to match with user's unit. The first one is au, the second is esu and last one is SI.311g
    for (auto it = alpha[frequency].begin(); it != alpha[frequency].end(); ++it) 
    {   
        if (unit == "au")
        {   
            replace(it->second[0].begin(), it->second[0].end(), 'D', 'E');
            temp.insert(make_pair(it->first, stod(it->second[0])));
        }
        else if (unit == "esu")
        {   
            replace(it->second[1].begin(), it->second[1].end(), 'D', 'E');
            temp.insert(make_pair(it->first, stod(it->second[1])));
        }
        else if (unit == "SI")
        {   
            replace(it->second[2].begin(), it->second[2].end(), 'D', 'E');
            temp.insert(make_pair(it->first, stod(it->second[2])));
        }
        else
        {
            throw runtime_error("ERROR in G16LOGfile::getAlpha(): Invalid unit. Please, use 'au', 'esu' or 'SI'.");
        }
        
    };

    return temp;
};

void G16LOGfile::setBeta()
{   
    string WhatPolar = "";
    vector<string> UsePolar;
    //auto UsePolar = this->vecPolarDip;
    for (int o = 0; o < 2; o++){
        if (o==0){
            WhatPolar = "input";
            UsePolar = this->vecPolarInp;
        }
        else{
            WhatPolar = "dipole";
            UsePolar = this->vecPolarDip;
        }
        map<double, int> start, end, start2, end2;
        for (int f = 0; f < this->vecNLOFrec.size(); f++)
        {
            double freqPrincipal = this->vecNLOFrec[f];
            for (int i = 0; i < UsePolar.size(); i++)
            {
                if (UsePolar[i].find("Beta(0;0,0)") != string::npos)
                {   
                    start.insert(make_pair(this->vecNLOFrec[0],i+2));
                    end.insert(make_pair(this->vecNLOFrec[0],i+18));
                    start2.insert(make_pair(this->vecNLOFrec[0],i+2));
                    end2.insert(make_pair(this->vecNLOFrec[0],i+18));
                }
                else if (UsePolar[i].find("Beta(-w;w,0) w= ") != string::npos)
                { 
                    string Freq = UsePolar[i].substr(UsePolar[i].find("Beta(-w;w,0)") + 15);
                    Freq = Freq.substr(0, Freq.find("nm"));
                    double FreqDouble = stod(Freq);
                    if(freqPrincipal == FreqDouble)
                    {
                        start.insert(make_pair(FreqDouble,i+2));
                        end.insert(make_pair(FreqDouble,i+26));
                    }
                }
                else if (UsePolar[i].find("Beta(-2w;w,w) w= ") != string::npos)
                { 
                    string Freq = UsePolar[i].substr(UsePolar[i].find("Beta(-2w;w,w)") + 16);
                    Freq = Freq.substr(0, Freq.find("nm"));
                    //convert Freq to double
                    double FreqDouble = stod(Freq);
                    if(freqPrincipal == FreqDouble)
                    {
                        start2.insert(make_pair(FreqDouble,i+2));
                        end2.insert(make_pair(FreqDouble,i+26));
                    }
                }

            }
        }
        map<double,map<string,vector<string>>> FrequencyInfo, FrequencyInfo2;
        map<string,vector<string>> NLOInfo;
        vector<vector<string>> Beta_0, Beta_w, Beta_2w;
        for (auto it = start.begin(); it != start.end(); ++it) 
        {   
            for(int i = start[it->first]; i < end[it->first]; i++)
            {   vector<string> T = this->customSplit(UsePolar[i]);
                if( T[0] == "||" && T[1] == "(z)"){
                    NLOInfo.insert(make_pair(T[0]+T[1], vector<string>({T[2], T[3], T[4]})));
                }   
                else{
                    NLOInfo.insert(make_pair(T[0], vector<string>({T[1], T[2], T[3]})));
                }         
            };
            FrequencyInfo.insert(make_pair(it->first, NLOInfo));
            NLOInfo.clear();
        }    
        this->Beta.insert(make_pair(WhatPolar, FrequencyInfo));

        for (auto it = start2.begin(); it != start2.end(); ++it) 
        {   
            for(int i = start2[it->first]; i < end2[it->first]; i++)
            {   vector<string> T = this->customSplit(UsePolar[i]);         
                if( T[0] == "||" && T[1] == "(z)"){
                    NLOInfo.insert(make_pair(T[0]+T[1], vector<string>({T[2], T[3], T[4]})));
                }   
                else{
                    NLOInfo.insert(make_pair(T[0], vector<string>({T[1], T[2], T[3]})));
                }         
            };
            FrequencyInfo2.insert(make_pair(it->first, NLOInfo));
            NLOInfo.clear();
        }    
        this->Beta2.insert(make_pair(WhatPolar, FrequencyInfo2));

    };
};


void G16LOGfile::setVibFrequencies()
{
    if ( this->vibFrequenciesStorage.size() == 0 )
    {
        this->vibFrequencies.emplace_back(0.0);
    };

    if ( this->vibFrequenciesStorage.size() > 0 )
    {
        for (int i = 0; i < this->vibFrequenciesStorage.size(); i++)
        {
            stringstream ss(this->vibFrequenciesStorage[i]);
            string line;
            while (getline(ss, line))
            {
                if (line.find("Frequencies --") != string::npos)
                {
                    istringstream iss(line);
                    vector<string> results((istream_iterator<string>(iss)), istream_iterator<string>());
                    for (int j = 2; j < results.size(); j++)
                    {
                        this->vibFrequencies.emplace_back(stod(results[j])); 
                    };
                };
            };
        };
    };
};

vector<double> G16LOGfile::getVibFrequencies()
{
    // cm^-1
    return this->vibFrequencies; // cm^-1
};

void G16LOGfile::setPrincipalAxesInertia()
{
    stringstream ss(this->principalAxisSTR);
    string line;

    while (getline(ss, line))
    {
        if (line.find("Eigenvalues --") != string::npos)
        {
            vector<string> results;
            {
                istringstream iss(line);
                for (string w; iss >> w; )
                    results.push_back(w);
            }

            if (results.size() == 5)
            {
                double ix = stod(results[2]);
                double iy = stod(results[3]);
                double iz = stod(results[4]);

                this->vecPrincipalAxesInertia.push_back(ix);
                this->vecPrincipalAxesInertia.push_back(iy);
                this->vecPrincipalAxesInertia.push_back(iz);
            }

            else
            {
                string raw = results[2];

                vector<string> base;
                {
                    size_t start = 0;
                    size_t pos;
                    while ((pos = raw.find('.', start)) != string::npos)
                    {
                        base.push_back(raw.substr(start, pos - start));
                        start = pos + 1;
                    }
                    base.push_back(raw.substr(start));
                }

                string inertia_x = base[0] + "." + base[1].substr(0, 5);
                string inertia_y = base[1].substr(5, 4) + "." + base[2].substr(0, 5);
                string inertia_z = base[2].substr(5, 4) + "." + base[3].substr(0, 5);

                double ix = stod(inertia_x); 
                double iy = stod(inertia_y); 
                double iz = stod(inertia_z); 

                this->vecPrincipalAxesInertia.push_back(ix);
                this->vecPrincipalAxesInertia.push_back(iy);
                this->vecPrincipalAxesInertia.push_back(iz);
            };
        };
    };
};

double G16LOGfile::getZPE()
{   
    // hartree
    return this->ZPE;
};

double G16LOGfile::getZPVE()
{   
    // kcal/mol
    return this->ZPVE;
};

double G16LOGfile::getH()
{
    // hartree
    return this->H;
};

double G16LOGfile::getG()
{
    // hartree
    return this->G;
};

double G16LOGfile::getS()
{
    // cal/(mol*K)
    return this->S;
};

double G16LOGfile::get_qEle()
{   
    return this->qEle;
};

void G16LOGfile::set_qVib()
{   
    double thetha_v = 0.0;
    double upper, lower;

    this->qVib = 1.0;

    if (this->vibFrequencies.size() == 0)
    {
        this->qVib = 1.0;
    }

    else
    {
        for (int i = 0; i < this->vibFrequencies.size(); i++)
        {
            if (this->vibFrequencies[i] <= 0)
            {
                continue;
            };
        
            thetha_v = (this->h * this->vibFrequencies[i] * this->c)/(this->pt.getConstant("k_B"));

            double upper = 1.0;
            double lower = 1.0 - exp(-thetha_v/(this->temperature));

            thetha_v = 0.0;
            this->qVib *= upper/lower;
        };
    };
};

void G16LOGfile::set_thetha_r()
{
    if (this->vecPrincipalAxesInertia.size() == 0)
    {
        this->thetha_r = 1.0;
    }

    else
    {
        double inertia_x = this->vecPrincipalAxesInertia[0]*this->pt.getConversion("amu_to_kg_m2"); // kg*m^2
        double inertia_y = this->vecPrincipalAxesInertia[1]*this->pt.getConversion("amu_to_kg_m2"); // kg*m^2
        double inertia_z = this->vecPrincipalAxesInertia[2]*this->pt.getConversion("amu_to_kg_m2"); // kg*m^2

        if (not this->isLinear)
        {
            double thetha_x = pow(this->h, 2) / (8 * pow(this->PI, 2) * this->pt.getConstant("k_B") * inertia_x);
            double thetha_y = pow(this->h, 2) / (8 * pow(this->PI, 2) * this->pt.getConstant("k_B") * inertia_y);
            double thetha_z = pow(this->h, 2) / (8 * pow(this->PI, 2) * this->pt.getConstant("k_B") * inertia_z);

            this->vecThetha_r = {thetha_x, thetha_y, thetha_z};
        }

        else
        {
            this->thetha_r = pow(this->h, 2) / (8 * pow(this->PI, 2) * this->pt.getConstant("k_B") * inertia_y);
        };
    };
};

void G16LOGfile::set_qRot()
{

    if (this->vecPrincipalAxesInertia.size() == 0)
    {
        this->qRot = 1.0;
    }

    else
    {
        if (this->isLinear)
        {
            {
                this->qRot = (this->temperature)/(this->sigma_r * this->thetha_r);
            }
        }
        else
        {
            double f1 = pow(this->PI, 0.5)/this->sigma_r;
            double f2 = pow(this->temperature, 1.5);
            double f3 = pow((this->vecThetha_r[0]*this->vecThetha_r[1]*this->vecThetha_r[2]), 0.5);

            this->qRot = f1 * (f2/f3);
        }
    };
};

void G16LOGfile::set_qTrans()
{
    double mass = this->mol.getMolecularMass()*this->pt.getConversion("amu_to_kg"); // kg

    double f1 = pow((2 * this->pt.getConstant("PI") * mass * this->pt.getConstant("k_B") * this->temperature) / (pow(this->pt.getConstant("h"), 2)), 1.5);
    double f2 = (this->pt.getConstant("k_B") * this->temperature) / this->pt.getConstant("P_0");

    this->qTrans = f1 * f2;
};

void G16LOGfile::set_qTot()
{
    this->qTot = this->qEle * this->qTrans * this->qRot * this->qVib;
};

double G16LOGfile::get_qTrans()
{
    return this->qTrans;
};

double G16LOGfile::get_qTot()
{
    return this->qTot;
};

double G16LOGfile::get_qRot()
{
    return this->qRot;
};

double G16LOGfile::get_qVib()
{
    return this->qVib;
};

void G16LOGfile::setIsLinear()
{
    int count = 0;
    for (int i = 0; i < this->vibFrequencies.size(); i++)
        {
            if (this->vibFrequencies[i] > 0)
            {
                count++;
            }
        }
    
    if (count > 1)
    {
        this->isLinear = false;
    }
    else
    {
        this->isLinear = true;
    }
};

map<string,double> G16LOGfile::getBeta(string orientation, string unit, double frequency, bool BSHG)
{   
    
    map<double,map<string,vector<string>>> beta;
    std::transform(orientation.begin(), orientation.end(), orientation.begin(),[](unsigned char c){ return std::tolower(c); });

    if (BSHG){
       beta = this->Beta2[orientation];
    }
    else
    {
        beta = this->Beta[orientation];
    }
    map<string,double> temp;
    if (!this->polarAsw)
    {
        throw runtime_error("ERROR in G16LOGfile::getBeta(): No NLO found in the log file., Try using polarAsw=1.");
    }
    //check if frequency is in the vecNLOFrec vector
    if (find(this->vecNLOFrec.begin(), this->vecNLOFrec.end(), frequency) == this->vecNLOFrec.end())
    {   
        string temp = "";
        for (int i = 0; i < this->vecNLOFrec.size(); i++)
        {
            temp += to_string(this->vecNLOFrec[i]) + ", ";
        }
        throw runtime_error("ERROR in G16LOGfile::getBeta(): Frequency not found in the log file. Try: " + temp + "instead.");
    }
    for (auto it = beta[frequency].begin(); it != beta[frequency].end(); ++it) 
    {   
        if (unit == "au")
        {   
            replace(it->second[0].begin(), it->second[0].end(), 'D', 'E');
            temp.insert(make_pair(it->first, stod(it->second[0])));
        }
        else if (unit == "esu")
        {   
            replace(it->second[1].begin(), it->second[1].end(), 'D', 'E');
            temp.insert(make_pair(it->first, stod(it->second[1])));
        }
        else if (unit == "SI")
        {   
            replace(it->second[2].begin(), it->second[2].end(), 'D', 'E');
            temp.insert(make_pair(it->first, stod(it->second[2])));
        }
        else
        {
            throw runtime_error("ERROR in G16LOGfile::getBeta(): Invalid unit. Please, use 'au', 'esu' or 'SI'.");
        }   
    };
    return temp;
};

void G16LOGfile::setGamma()
{   
    string WhatPolar = "";
    vector<string> UsePolar;
    //auto UsePolar = this->vecPolarDip;
    for (int o = 0; o < 2; o++){
        if (o==0){
            WhatPolar = "input";
            UsePolar = this->vecPolarInp;
        }
        else{
            WhatPolar = "dipole";
            UsePolar = this->vecPolarDip;
        }
        map<double, int> start, end, start2, end2;
        for (int f = 0; f < this->vecNLOFrec.size(); f++)
        {
            double freqPrincipal = this->vecNLOFrec[f];
            for (int i = 0; i < UsePolar.size(); i++)
            {
                if (UsePolar[i].find("Gamma(0;0,0,0)") != string::npos)
                {   
                    
                    start.insert(make_pair(this->vecNLOFrec[0],i+2));
                    end.insert(make_pair(this->vecNLOFrec[0],i+19));
                    start2.insert(make_pair(this->vecNLOFrec[0],i+2));
                    end2.insert(make_pair(this->vecNLOFrec[0],i+19));
                }
                else if (UsePolar[i].find("Gamma(-w;w,0,0) w= ") != string::npos)
                { 
                    string Freq = UsePolar[i].substr(UsePolar[i].find("Gamma(-w;w,0,0)") + 18);
                    Freq = Freq.substr(0, Freq.find("nm"));
                    double FreqDouble = stod(Freq);
                    if(freqPrincipal == FreqDouble)
                    {
                        start.insert(make_pair(FreqDouble,i+2));
                        end.insert(make_pair(FreqDouble,i+40));
                    }
                }
                else if (UsePolar[i].find("Gamma(-2w;w,w,0) w= ") != string::npos)
                { 
                    string Freq = UsePolar[i].substr(UsePolar[i].find("Gamma(-2w;w,w,0)") + 19);
                    Freq = Freq.substr(0, Freq.find("nm"));
                    //convert Freq to double
                    double FreqDouble = stod(Freq);
                    if(freqPrincipal == FreqDouble)
                    {
                        start2.insert(make_pair(FreqDouble,i+2));
                        end2.insert(make_pair(FreqDouble,i+58));
                    }
                }

            }
        }
        map<double,map<string,vector<string>>> FrequencyInfo, FrequencyInfo2;
        map<string,vector<string>> NLOInfo;
        //vector<vector<string>> Gamma_0, Gamma_w, Gamma_2w;
        for (auto it = start.begin(); it != start.end(); ++it) 
        {   
            for(int i = start[it->first]; i < end[it->first]; i++)
            {   vector<string> T = this->customSplit(UsePolar[i]);
                NLOInfo.insert(make_pair(T[0], vector<string>({T[1], T[2], T[3]})));        
            };
            FrequencyInfo.insert(make_pair(it->first, NLOInfo));
            NLOInfo.clear();
        }    
        this->Gamma.insert(make_pair(WhatPolar, FrequencyInfo));

        for (auto it = start2.begin(); it != start2.end(); ++it) 
        {   
            for(int i = start2[it->first]; i < end2[it->first]; i++)
            {   vector<string> T = this->customSplit(UsePolar[i]);         
                NLOInfo.insert(make_pair(T[0], vector<string>({T[1], T[2], T[3]})));       
            };
            FrequencyInfo2.insert(make_pair(it->first, NLOInfo));
            NLOInfo.clear();
        }    
        this->Gamma2.insert(make_pair(WhatPolar, FrequencyInfo2));

   };
// Print this->Gamma2 on terminal
};

map<string,double> G16LOGfile::getGamma(string orientation, string unit, double frequency, bool GSHG)
{   
    std::transform(orientation.begin(), orientation.end(), orientation.begin(),[](unsigned char c){ return std::tolower(c); });

    map<double,map<string,vector<string>>> gamma;
    if (GSHG){
       gamma = this->Gamma2[orientation];
    }
    else
    {
        gamma = this->Gamma[orientation];
    }
    map<string,double> temp;
    if (!this->polarAsw)
    {
        throw runtime_error("ERROR in G16LOGfile::getGamma(): No NLO found in the log file., Try using polarAsw=1.");
    }
    //check if frequency is in the vecNLOFrec vector
    if (find(this->vecNLOFrec.begin(), this->vecNLOFrec.end(), frequency) == this->vecNLOFrec.end())
    {   
        string temp = "";
        for (int i = 0; i < this->vecNLOFrec.size(); i++)
        {
            temp += to_string(this->vecNLOFrec[i]) + ", ";
        }
        throw runtime_error("ERROR in G16LOGfile::getGamma(): Frequency not found in the log file. Try: " + temp + "instead.");
    }
    for (auto it = gamma[frequency].begin(); it != gamma[frequency].end(); ++it) 
    {   
        if (unit == "au")
        {   
            replace(it->second[0].begin(), it->second[0].end(), 'D', 'E');
            temp.insert(make_pair(it->first, stod(it->second[0])));
        }
        else if (unit == "esu")
        {   
            replace(it->second[1].begin(), it->second[1].end(), 'D', 'E');
            temp.insert(make_pair(it->first, stod(it->second[1])));
        }
        else if (unit == "SI")
        {   
            replace(it->second[2].begin(), it->second[2].end(), 'D', 'E');
            temp.insert(make_pair(it->first, stod(it->second[2])));
        }
        else
        {
            throw runtime_error("ERROR in G16LOGfile::getGamma(): Invalid unit. Please, use 'au', 'esu' or 'SI'.");
        }   
    };
    return temp;
};


vector<string> G16LOGfile::customSplit(string str, char separator) 
{   
    vector<string> tokens;
    istringstream iss(str);
    string token;

    while (getline(iss, token, separator)) 
    {
        if (!token.empty()) {  // Skip empty tokens
            tokens.push_back(token);
        }
    };
    return tokens;
};

string G16LOGfile::toStr()
{
    return "G16LOGFile: Calculation of " + this->str_filePath.substr(this->str_filePath.find_last_of("/") + 1) + " done in " + this->info + ", with the level of theory " + this->method + "/" + this->basisValue + " and " + "SCF energy of " + to_string(this->scfValue) + " Hartrees.";
};




//!----------------------- Destructor -----------------------//

G16LOGfile::~G16LOGfile()
{
    //clear int
    this->charge = 0;
    this->multiplicity = 0;

    // clear the size_t
    this->dipoleFinder = 0;
    this->scf = 0;
    this->scfC = 0;
    this->aHomoFinder = 0;
    this->aLumoFinder = 0;
    this->tdFinder = 0;
    this->polarFinder = 0;
    this->normalT = 0;
    this->stdT = 0;
    this->starterSCF = 0;
    this->endSCF = 0;
    this->starterMethod = 0;
    this->basis = 0;

    // clear the strings
    this->moleculeSTR.clear();
    //this->polarAuxiliary.clear();

    // clear the doubless
    this->scfValue = 0;
    this->aHomoValue = 0;
    this->aLumoValue = 0;
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

//!----------------------- Notepad -----------------------//

//! TODO:

// TODO:getAlpha
// TODO:getBeta
// TODO:getGamma
// TODO:getGradient
// TODO:getTransitions_str    // Create a vector to store the lines of the stringstream

// TODO:getWavelengths
// TODO:getSymmetries
// TODO:getSymmetry
// TODO:getOscillatorForce
// TODO:getTransContributions
// TODO:destructor for streamstrings sdip and sinp

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


//! BUGS:

// Wrong HOMO for NLO calculations; Need to make a storage like Molecule;
