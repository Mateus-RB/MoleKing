//   MoleKing //
//
//   File:        [G16LOGfile.hpp]
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

#ifndef G16Process_hpp
#define G16Process_hpp

#include <stdio.h>
#include <string>
#include <fstream>
#include "../myMath/Vectors.hpp"
#include "../myMath/Matrix.hpp"
#include <iostream>
#include <numeric>
#include <regex>
#include "../chemicalUnits/Molecule.hpp"
#include <stdexcept>
#include <map>
#include <cmath>

using namespace std;
#endif /* Testing_hpp */

class G16LOGfile{
private:

    //* Constants

    const double k_B = 1.3806503E-23;                    // Boltzmann (J/K)
    const double h = 6.62607004E-34;                     // Planck (J*s)
    const double hbar = 1.0545718e-34;                   // Planck reduced (J*s)
    const double c = 2.99792458e10;                      // Speed of light (cm/s)
    const double R = 1.98720425864083;                   // Gas constant (cal/(mol*K))
    const double PI = 3.14159265358979323846;            // Pi
    const double P_0  = 101325;                          // Pressure (Pa)
    const double AMU_2_KG = 1.660538921E-27;             // Conversion from AMU to Kg
    const double AMU_2_KG_M2 = 4.65082513926e-48;        // Conversion from AMU to Kg/m^2  


    //* istringstream

    istringstream logfile;

    //* int

    int charge;
    int multiplicity;
    int link;
    int countNormal;


    //* size_t
    size_t dipoleFinder;
    size_t mullikenFinder;
    size_t aHomoFinder;
    size_t aLumoFinder;
    size_t bHomoFinder;
    size_t bLumoFinder;
    size_t tdFinder;
    size_t scf;
    size_t scfC;
    size_t normalT;   
    size_t stdT;
    size_t cpFinder;
    size_t starterSCF;
    size_t endSCF;
    size_t starterMethod;  
    size_t basis; 
    size_t chargeMultiFinder;
    size_t polarFinder;

    //* string
    string basisValue;    
    string line;
    string value;
    string moleculeRange;
    string info;
    string method;
    string moleculeSTR = "";
    string cpSTR = "";
    string polarSTR = "";
    string vibFreqSTR = "";
    string mullikenSTR = "";
    string str_filePath;
    string linkStorageSTD = "";
    string polarAuxiliaryDip;
    string polarAuxiliaryInp;
    string Freq;
    string dof;
    string principalAxisSTR;
    vector <string> iptStorage;
    vector <string> stdStorage;
    vector <string> aHomoStorage;
    vector <string> aLumoStorage;
    vector <string> bHomoStorage;
    vector <string> bLumoStorage;
    vector <string> dipoleStorage;
    vector <string> tdStorage;
    vector <string> polarStorageDip;
    vector <string> polarStorageInp;
    vector <string> vecPolarDip;
    vector <double> vibFrequencies;
    vector <string> vecPolarInp;
    vector <string> linkStorage;
    vector <string> vibFrequenciesStorage;
    vector<string> elstDipoleStorage;
    vector<string> alphaStorage;
    vector<string> betaStorage;
    vector<string> gammaStorage;

    //* double
    double temperature;
    double scfValue;
    double ZPE;
    double ZPVE;
    double aHomoValue;
    double aLumoValue;
    double bHomoValue;
    double bLumoValue;
    double dipoleTot;
    double dipoleX;
    double dipoleY;
    double dipoleZ;
    double FreqDouble;   
    double sigma_r;
    double H;
    double S;
    double G;
    double qEle;
    double qVib;
    double qRot;
    double thetha_r;
    double qTrans;
    double qTot;

    //* molecule
    Molecule mol;

    //* bool

    bool ntFound;
    bool stdFound;
    bool scfConvergence;
    bool polarAsw;
    bool tdAsw;
    bool thermoAsw;
    bool cpAsw;
    bool dipFinder = false;
    bool isLinear = false;

    //* map and vectors

    map<string, vector<string>> Orbitals;
    map<string, vector<string>> bOrbitals;
    map<int, map<string, double>> transitions;
    map<string,map<double,map<string,vector<string>>>> Alpha;
    map<string,map<double,map<string,vector<string>>>> Beta;
    map<string,map<double,map<string,vector<string>>>> Beta2;
    map<string,map<double,map<string,vector<string>>>> Gamma;
    map<string,map<double,map<string,vector<string>>>> Gamma2;

    vector<string> Occupied;
    vector<string> bOccupied;
    vector<string> Unoccupied;
    vector<string> bUnoccupied;
    vector<double> vecNLOFrec;
    vector<double> vecFrecDip;
    vector<double> vecFrecInp;
    vector<double> vecPrincipalAxesInertia;
    vector<double> vecThetha_r;
    
    //* teste

    string aLumoStorageSTR = "";
    string aHomoStorageSTR = "";
    string bLumoStorageSTR = "";
    string bHomoStorageSTR = "";

    //* set functions
    void detectLink();
    void readLOGFile();
    void setMolecule();
    void setChargePoints();
    void setOrbitals();
    void setHOMO();
    void setLUMO();
    void setTransitions();
    void setDipole();
    void splitter();
    void setNLO();
    void setNLOFrequency();
    void setVibFrequencies();
    void setAlpha();
    void setBeta();
    void setIsLinear();
    void setSigma_r();
    void setPrincipalAxesInertia();
    void setGamma();
    void set_qVib();
    void set_thetha_r();
    void set_qRot();
    void set_qTrans();
    void set_qEle();
    void set_qTot();
    vector<string> customSplit(string str, char separator = ' ');

public:
    G16LOGfile(string filePath, bool polarAsw = 0, bool tdAsw = 0, bool thermoAsw = 0, bool cpAsw = 0, double temperature = 298.15, int link = -1);   
    ~G16LOGfile();      
    map<string, vector<string>> getOrbitals(); 
    map<int, map<string, double>> getTransitions(const int index = 0);
    double getEnergy();
    double getZPE();
    double getZPVE();
    double getH();
    double getS();
    double getG();
    double get_qEle();
    double get_qVib();
    double get_qRot();
    double get_qTrans();
    double get_qTot();
    vector<double> getHOMO(int index = -1);
    vector<double> getLUMO(int index = 0);    
    double getDipole(string axis = "tot");
    string getDate();    
    string getBasis();
    string getMethod();
    map<string,double> getAlpha(string orientation = "Dipole", string unit = "esu", double frequency = 0);
    map<string,double> getBeta(string orientation = "Dipole", string unit = "esu", double frequency = 0, bool BSHG = 0);
    map<string,double> getGamma(string orientation = "Dipole", string unit = "esu", double frequency = 0, bool GSHG = 0);
    vector<double> getVibFrequencies();
    string toStr();
    Molecule getMolecule();
    vector<double> getNLOFrequency();  
    vector<string> getNLO(string orientation = "input");
    
};
