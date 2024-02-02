//   CPP FILE HEADER //
//
//   File:        [G16Process.cpp]
//
//   Author(s):   ['LEEDMOL Reserch Group']
//   Site(s):     ['https://www.researchgate.net/lab/LEEDMOL-Heibbe-Cristhians-Research-Group-Heibbe-Cristhian-B-De-Oliveira']
//   Email(s):    ['heibbe@ufg.br']
//   Credits:     ['Copyright © 2023 LEEDMOL. All rights reserved.']
//   Date:        ['17.01.2023']
//   Version:     ['1.5.0']
//   Status:      ['Development']
//   Language:    ['C++','Python']
//   Description: ['A python module written in C++ for theoretical chemestry']

#include "G16Process.hpp"

//!----------------------- G16LOGfile -----------------------//

G16LOGfile::G16LOGfile(string filePath, bool polarAsw, bool tdAsw, int link)
    : ntFound(false), stdFound(false), scfConvergence(true),str_filePath(filePath), polarAsw(polarAsw), tdAsw(tdAsw)
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
            cerr << "WARNING in G16LOGfile::G16LOGfile(): We detected " << this->linkStorage.size() << " link calculation(s) in your log file. As default, we will use the last one. To change it use: G16LOGfile(filename, link=n)" << endl;
            this->logfile.str(this->linkStorage[0]);
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
        this->vecFrec = {0.000};
        setNLO();
        setFrequency();
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
        dipoleFinder = line.find("Tot=");
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

            this->scfValue = stod(value.substr(starterSCF + 3, endSCF - starterSCF - 3));
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
                    aLumoStorageSTR += line + "\n";
                    break;
                };
                aHomoStorageSTR += line + "\n";
            };
            this->aHomoStorage.emplace_back(aHomoStorageSTR);
            aHomoStorageSTR = "";
        };
        if (aLumoFinder != string::npos)
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
        if (line.find("Input orientation:") != string::npos)
        {
            while (getline(this->logfile, line))
            {
                if (line.find("Distance matrix (angstroms):") != string::npos)
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
            moleculeSTR = moleculeSTR.substr(0, moleculeSTR.rfind("\n"));
            iptStorage.emplace_back(moleculeSTR);
            moleculeSTR = "";
        };
        // If the line contains "Standard orientation:", extract the molecule's geometry in Standard Orientation
        if (line.find("Standard orientation:") != string::npos)
        {
            while (getline(this->logfile, line))
            {
                if (line.find(" Rotational constants (GHZ): ") != string::npos)
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
            moleculeSTR = moleculeSTR.substr(0, moleculeSTR.rfind("\n"));
            stdStorage.emplace_back(moleculeSTR);
            moleculeSTR = "";
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
        if (polarAsw)
        {
            if (line.find(" Dipole orientation:") != string::npos)
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
                polarStorage.emplace_back(polarSTR);
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
        cerr << "WARNING in G16LOGfile::getEnergy(): SCF convergence not achieved. Please check your log file for 'Convergence criterion not met.'" << endl;
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
        //cerr << "WARNING in G16LOGfile::getMolecule(): The geometry of " << temp << " was not taken from the standard orientation." << endl;
    };
    return this->mol;
};

// Function to user get the orbitals
map<string, vector<string>> G16LOGfile::getOrbitals() 
{   
    if (this->bOccupied.size() > 0 && this->bUnoccupied.size() > 0)
    {   
        //throw an warning
        cerr << "WARNING in G16LOGfile::getOrbitals(): Your job is an open shell. Please, use getOrbitals()['Alpha_Occupied'] or ['Beta_Occupied']" << endl;
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
    // get the last element of the polarStorage vector
    this->polarAuxiliary = this->polarStorage[this->polarStorage.size() - 1];
    // create a stringstream object with polarAuxiliary
    stringstream ss(polarAuxiliary);
    // Store each line of the stringstream in the lines vector
    while (getline(ss, line))
    {
        this->vecPolar.emplace_back(line);
    };
    //check if the last element of the vector lines is equal to "" or " " or "\n"
    if (this->vecPolar[vecPolar.size() - 1] == "" || this->vecPolar[vecPolar.size() - 1] == " " || this->vecPolar[vecPolar.size() - 1] == "\n")
    {
        // if true, remove the last element of the vector lines
        this->vecPolar.pop_back();
    };

    
};


void G16LOGfile::setFrequency()
{    
    for (int i = 0; i < this->vecPolar.size(); i++)
    {
        if (this->vecPolar[i].find("Alpha(-w;w)") != string::npos)
        {
            this->Freq = this->vecPolar[i].substr(this->vecPolar[i].find("Alpha(-w;w)") + 16);
            this->Freq = this->Freq.substr(0, this->Freq.find("nm"));
            //convert Freq to double
            this->FreqDouble = stod(this->Freq);
            // add the frequency to the vecFrec vector
            this->vecFrec.emplace_back(this->FreqDouble);
        };
    };
};

vector<double> G16LOGfile::getFrequency()
{   
    if (!this->polarAsw)
    {
        throw runtime_error("ERROR in G16LOGfile::getFrequency(): No frequency found in the log file., Try using polarAsw=1.");
    }
    return this->vecFrec;
    
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
