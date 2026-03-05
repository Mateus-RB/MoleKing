//   CPP FILE HEADER //
//
//   File:        [OrcaProcess.cpp]
//
//   Author(s):   ['LEEDMOL Reserch Group']
//   Site(s):     ['https://www.researchgate.net/lab/LEEDMOL-Heibbe-Cristhians-Research-Group-Heibbe-Cristhian-B-De-Oliveira']
//   Email(s):    ['heibbe@ufg.br']
//   Credits:     ['Copyright © 2023 LEEDMOL. All rights reserved.']
//   Date:        ['05.04.2026']
//   Version:     ['1.6.0']
//   Status:      ['Development']
//   Language:    ['C++','Python']
//   Description: ['A python module written in C++ for theoretical chemistry']

#include "OrcaProcess.hpp"

//!----------------------- OrcaLOGfile -----------------------//

ORCALOGfile::ORCALOGfile(string filePath, bool thermoAsw)
    : str_filePath(filePath), thermoAsw(thermoAsw)
{

    // SET FUNCTIONS
    
    initializeLOGFile();
    readLOGFile();
    setMolecule();

    if (thermoAsw)
    {
        setVibFrequencies();
        setIsLinear();
    //     setPrincipalAxesInertia();
    //     // set_qVib();
    //     // set_qRot();
    //     // set_qTrans();
    //     // set_qTot();
    };
}

void ORCALOGfile::initializeLOGFile()
{
    // Read file and check if it exists
    ifstream log_file(this->str_filePath);

    if (!log_file.is_open())
    {
        throw runtime_error("ERROR in ORCALOGfile::ORCALOGfile(): File not found. Please check your file path.");
    }

    // Read entire file content into a string
    string fileContent;
    string lineBuffer;
    while (getline(log_file, lineBuffer))
    {
        fileContent += lineBuffer + "\n";
    }
    log_file.close();

    // Set istringstream with file content
    this->logfile.str(fileContent);

}

void ORCALOGfile::readLOGFile()
{
    while (getline(this->logfile, line))
    {   
        if (line.find("CARTESIAN COORDINATES (ANGSTROEM)") != string::npos)
        {
            moleculeSTR = "";  // Clear previous content
            while (getline(this->logfile, line))
            {
                if (line.find("CARTESIAN COORDINATES (A.U.)") != string::npos)
                {
                    break;
                };
                moleculeSTR += line + "\n";
            }
            // removing the first line of the moleculeSTR, which is the header of the geometry block. and also removing the last two lines, which are empty lines.
            for (int i = 0; i < 1; i++) 
            {
                moleculeSTR = moleculeSTR.substr(moleculeSTR.find("\n") + 1);
            };
            moleculeSTR = moleculeSTR.substr(0, moleculeSTR.rfind("----------------------------"));
            moleculeSTR = moleculeSTR.substr(0, moleculeSTR.rfind("\n"));
            
            geometryBlocks.emplace_back(moleculeSTR);
            moleculeSTR = ""; // Clear moleculeSTR for the next geometry block
        }
        if (line.find(" Total Charge") != string::npos)
        {   
            //  Total Charge           Charge          ....    0
            istringstream iss(line);
            vector<string> tokens((istream_iterator<string>(iss)), istream_iterator<string>());
            this->charge = stoi(tokens.back());
        }
        if (line.find("Multiplicity") != string::npos)
        {   
            istringstream iss(line);
            vector<string> tokens((istream_iterator<string>(iss)), istream_iterator<string>());
            this->multiplicity = stoi(tokens.back());
            this->qEle = stoi(tokens.back());
        }
        
        if (line.find("FINAL SINGLE POINT ENERGY") != string::npos)
        {   
            istringstream iss(line);
            vector<string> tokens((istream_iterator<string>(iss)), istream_iterator<string>());
            this->scfValue = stod(tokens.back()); // Hartree
        }

        if (thermoAsw)
        {
            if (line.find("Scaling factor for frequencies = ") != string::npos)
            {
                getline(this->logfile, line);
                while (getline(this->logfile, line))
                {
                    if (line.empty())
                    {
                        break;
                    };
                    this->vibFrequenciesStorage.emplace_back(line); // cm**-1
                }
            }
            if (line.find("Zero point energy") != string::npos)
            {
                istringstream iss(line);
                vector<string> tokens((istream_iterator<string>(iss)), istream_iterator<string>());
                // Zero point energy                ...      0.09468478 Eh      59.42 kcal/mol
                // The ZPE in hartree is the [4] element of the tokens vector.
                this->ZPE = stod(tokens[4]); // Hartree
            }
            if (line.find("Total enthalpy") != string::npos)
            {
                istringstream iss(line);
                vector<string> tokens((istream_iterator<string>(iss)), istream_iterator<string>());
                // Total enthalpy                    ...      0.09468478 Eh      59.42 kcal/mol
                // The enthalpy in hartree is the [3] element of the tokens vector.
                this->H = stod(tokens[3]); // Hartree
            }
            if (line.find("Final Gibbs free energy ") != string::npos)
            {
                istringstream iss(line);
                vector<string> tokens((istream_iterator<string>(iss)), istream_iterator<string>());
                // Final Gibbs free energy         ...   -152.06333403 Eh
                // The enthalpy in hartree is the [5] element of the tokens vector.
                this->G = stod(tokens[5]); // Hartree                
            }
            if (line.find("Final entropy term                ... ") != string::npos)
            {
                istringstream iss(line);
                vector<string> tokens((istream_iterator<string>(iss)), istream_iterator<string>());
                // Final entropy term                ...      0.03038983 Eh     19.07 kcal/mol
                // The enthalpy in hartree is the [6] element of the tokens vector.
                this->S = (stod(tokens[6])/298.15)*1000; // cal/mol.K                                
            }
            if (line.find("Symmetry Number:") != string::npos)
            {
                istringstream iss(line);
                vector<string> tokens((istream_iterator<string>(iss)), istream_iterator<string>());
                // Point Group:  C1, Symmetry Number:   1  
                // The sigma_r is the last element of the tokens vector.
                this->sigma_r = stod(tokens.back()); // dimensionless
            }
            if (line.find("Rotational constants in cm-1: ") != string::npos)
            {
                istringstream iss(line);
                vector<string> tokens((istream_iterator<string>(iss)), istream_iterator<string>());
                // Rotational constants in cm-1:             1.165006             0.304730             0.265422 
                // or
                // Rotational constants in cm-1:             0.000000             0.000000             0.000000 

                if (stod(tokens[4]) != 0.0)
                {
                    double temp_thethar_x = stod(tokens[4]) * this->pt.getConversion("cm1_to_K"); // Convert from cm^-1 to K
                    double temp_thethar_y = stod(tokens[5]) * this->pt.getConversion("cm1_to_K"); // Convert from cm^-1 to K
                    double temp_thethar_z = stod(tokens[6]) * this->pt.getConversion("cm1_to_K"); // Convert from cm^-1 to K
                    this->vecThetha_r.push_back(temp_thethar_x);
                    this->vecThetha_r.push_back(temp_thethar_y);
                    this->vecThetha_r.push_back(temp_thethar_z);
                }
                else
                {
                    this->thetha_r = 1.0; // If the rotational constant is zero, it means that the molecule is linear and the rotational partition function should be divided by the symmetry number only, without the sqrt(pi*thetha_r) term.
                }
            }
        }
    };

    // The molecule's geometry always is the last geometry!
    
    this->moleculeSTR = geometryBlocks[geometryBlocks.size() - 1];
};

void ORCALOGfile::setMolecule()
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
    // If the molecule object has no atoms, throw an error
    if (this->mol.getSize() == 0)
    {
        throw runtime_error("Molecule not found in the log. Please check your output file.");
    };
    this->mol.setCharge(this->charge);
    this->mol.setMultiplicity(this->multiplicity);
};

void ORCALOGfile::setIsLinear()
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

void ORCALOGfile::setVibFrequencies()
{
    if (this->vibFrequenciesStorage.size() == 0 )
    {
        this->vibFrequencies.emplace_back(0.0);
    };

    if ( this->vibFrequenciesStorage.size() > 0 )
    {
        for (int i = 0; i < this->vibFrequenciesStorage.size(); i++)
        {   
            istringstream iss(this->vibFrequenciesStorage[i]);
            vector<string> results((istream_iterator<string>(iss)), istream_iterator<string>());

            if (stod(results[1]) != 0.0)
            {
                this->vibFrequencies.emplace_back(stod(results[1]));       
            }     
        }
    }

    if (this->vibFrequencies.size() == 0)
    {
        this->vibFrequencies.emplace_back(0.0);
    };
};

void ORCALOGfile::set_qVib(double Temperature)
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
        
            thetha_v = (this->pt.getConstant("h") * this->vibFrequencies[i] * this->pt.getConstant("C"))/(this->pt.getConstant("k_B"));

            double upper = 1.0;
            double lower = 1.0 - exp(-thetha_v/(Temperature));

            thetha_v = 0.0;
            this->qVib *= upper/lower;
        };
    };
};

void ORCALOGfile::set_qRot(double Temperature)
{
    if (this->mol.getSize() == 1)
    {
        this->qRot = 1.0;
    }

    else
        {
        if (this->isLinear)
        {
            {   
                //this->qRot = (Temperature)/(this->thetha_r);
                this->qRot = (Temperature)/(this->sigma_r * this->thetha_r);
            }
        }
        else
        {
            double f1 = pow(this->pt.getConstant("PI"), 0.5)/this->sigma_r;
            //double f1 = sqrt(pow(this->PI, 0.5));
            double f2 = pow(Temperature, 1.5);
            double f3 = pow((this->vecThetha_r[0]*this->vecThetha_r[1]*this->vecThetha_r[2]), 0.5);

            this->qRot = f1 * (f2/f3);
        }
    }
};

void ORCALOGfile::set_qTrans(double Temperature)
{
    double mass = this->mol.getMolecularMass()*this->pt.getConversion("amu_to_kg"); // kg

    double f1 = pow((2 * this->pt.getConstant("PI") * mass * this->pt.getConstant("k_B") * Temperature) / (pow(this->pt.getConstant("h"), 2)), 1.5);
    double f2 = (this->pt.getConstant("k_B") * Temperature) / this->pt.getConstant("P_0");

    this->qTrans = f1 * f2;
};

double ORCALOGfile::get_qVib(double Temperature)
{
    this->set_qVib(Temperature);
    return this->qVib;
};

void ORCALOGfile::set_qTot()
{
    this->qTot = this->qEle * this->qTrans * this->qRot * this->qVib;
};


double ORCALOGfile::get_qRot(double Temperature)
{
    this->set_qRot(Temperature);
    return this->qRot;
};

double ORCALOGfile::get_qTrans(double Temperature)
{
    this->set_qTrans(Temperature);
    return this->qTrans;
};

double ORCALOGfile::get_qTot(double Temperature)
{
    this->get_qEle();
    this->set_qVib(Temperature);
    this->set_qRot(Temperature);
    this->set_qTrans(Temperature);
    this->set_qTot();
    return this->qTot;
};

vector<double> ORCALOGfile::getVibFrequencies()
{
    // cm^-1
    return this->vibFrequencies; // cm^-1
};

double ORCALOGfile::getEnergy()
{
    return this->scfValue; // Hartree
};

double ORCALOGfile::getZPE()
{
    return this->scfValue + this->ZPE; // Hartree
};

double ORCALOGfile::get_qEle()
{   
    return this->qEle; // dimensionless
};

double ORCALOGfile::get_sigmaR()
{
    return this->sigma_r; // Dimensioless
};

double ORCALOGfile::getH()
{
    return this->H; // Hartree  
};

double ORCALOGfile::getG()
{
    return this->G; // Hartree  
};

double ORCALOGfile::getS()
{
    return this->S; // Kcal/mol.K
};


Molecule ORCALOGfile::getMolecule()
{
    if (this->mol.getSize() == 0)
    {
        throw runtime_error("Molecule not found in the log. Please check your output file.");
    }
    return this->mol;
};

ORCALOGfile::~ORCALOGfile() noexcept
{
    this->str_filePath = "";
    this->line = "";
    this->moleculeSTR = "";
};