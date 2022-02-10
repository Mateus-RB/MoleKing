//
//  G16Process.hpp
//  MoleKing_util
//
//  Created by Thiago Lopes on 01/03/20.
//  Copyright © 2020 Laboratório de Estrutura Eletrônica e Dinâmica Molecular. All rights reserved.
//

#ifndef G16Process_hpp
#define G16Process_hpp

#include <stdio.h>
#include <string>
#include <fstream>
#include "../math/Vectors.hpp"
#include "../math/Matrix.hpp"
#include <iostream>
#include <regex>
#include "../chemicalUnits/Molecule.hpp"
#include <stdexcept>
#include <map>

using namespace std;
#endif /* G16Process_hpp */

class ExcStates{
private:
    int statesNumber;
    vector<double> wlValues;
    vector <double> energies;
    vector <double> oscillator;
    vector <string> symmetries;
    vector <vector <pair < pair <int, int>, double > > > transitions;
    
public:
    ExcStates(int statesNumber);
    ExcStates();
    
    void setWavelength(int state, double value);
    void setEnergy(int state, double value);
    void setOscillatorForce(int state, double value);
    void setTransitions(int state, vector <pair < pair <int, int>, double > > values);
    void setSymmetry(int state, string value);
    
    int getstatesNumber();
    string getSymmetry(int state);
    double getWavelength(int state);
    double getEnergy(int state);
    double getOscillatorForce(int state);
    vector <pair < pair <int, int>, double > > getTransition(int state);
    vector < pair <string, double> > getTransContribution(int state);
};

class PolarValues{
private:
    vector<string> dName;
    vector< pair <string, vector<string> > > aName, bName, gName;
    vector<double> dValue;
    vector< pair <string, vector<double> > > aValue, bValue, gValue;
    
public:
    PolarValues();
    
    void setDipole(string name, double value);
    void setAlpha(string eleName, string name, double value);
    void setBeta(string eleName, string name, double value);
    void setGamma(string eleName, string name, double value);
    
    double getDipole(string name);
    double getAlpha(string eleName, string name);
    double getBeta(string eleName, string name);
    double getGamma(string eleName, string name);
    
};

class GradientValues{
private:
    Matrix cartesianGradient;
public:
    GradientValues(); 

    void setGradient(Matrix gradient);
    Matrix getGradient();
};

class G16LOGfile{
private:
    double energy;
    string filePath, fileType, levelTheory, basis, date;
    int size;
    bool polarAsw, optAsw, stateAsw, calcDone, not_stop, chelpg;
    vector<double> occOrb, virtOrb;
    Molecule molecule;
    PolarValues polarValues;
    GradientValues gradientValues;
    void makePolar(vector <string> fileLines);
    void makeStates(vector <string> fileLines);
    void molConstructor(vector <string> fileLines);
    void makeGradient(vector <string> fileLines);
    int statesNum(vector <string> fileLines);
    ExcStates exSates;
    vector <string>  getTransition(int state);
    
public:
    G16LOGfile(string filePath, bool polarAsw = 0);
    
    double scfEnergy();
    Molecule getMolecule();
    double getDipole(string name);
    double getAlpha(string eleName, string name);
    double getBeta(string eleName, string name);
    double getGamma(string eleName, string name);
    string toStr();
    double getOscillatorForce(int state);
    double getWavelength(int state);
    string getSymmetry(int state);
    vector <double> getOscillatorForces();
    vector <double> getWavelengths();
    vector <string> getSymmetries();
    vector <vector <string> > getTransitions();
    vector <string> getTransitionsStr();
    vector <vector <string>> getTransContributions();
    Matrix getGradient();
};

class G16FCHKfile{
private:
    double energy;
    string filePath, fileType, levelTheory, basis, date;
    int size, electronNumber;
    bool  optAsw, calcDone;
    Matrix cartesianGradient, quadrupoleMoment;
    vector<double> occOrb, virtOrb;
    GradientValues gradientValues;
    Molecule molecule;
    void molConstructor(vector <string> fileLines);
    void makeGradient(vector <string> fileLines);

public:
    G16FCHKfile(string filePath);
    Molecule getMolecule();
    Matrix getCartesianGradient();
};

