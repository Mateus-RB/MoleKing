//
//  G16Process.cpp
//  MoleKing_util
//
//  Created by Thiago Lopes on 01/03/20.
//  Copyright © 2020 Laboratório de Estrutura Eletrônica e Dinâmica Molecular. All rights reserved.
//

#include "G16Process.hpp"

vector <string> splitString(string lineSTR, char splitTarget){
    vector <string> splittedLine;
    string temp;
    vector<char> lineVector(lineSTR.begin(),lineSTR.end());
    for (int i = 0; i < (int) lineVector.size(); i++){
        if (lineVector[i] != splitTarget){
            temp = temp + lineVector[i];
        } else {
            if (temp != ""){
                splittedLine.push_back(temp);
                temp = "";
            };
        };
    };
    if (temp != "" || temp != " "){
        splittedLine.push_back(temp);
    }
    return splittedLine;
};

/*
----------------------- G16LOGfile -----------------------
*/

G16LOGfile::G16LOGfile(string filePath, bool polarAsw){
    ifstream arq;
    this->polarAsw = polarAsw;
    arq.open(filePath, ifstream::in);
    string lineSTR;
    regex states_re("(.*)Excitation energies and oscillator strengths:");
    regex done_re(" Normal termination of (.*)");
    regex opt_re("(.*) opt(.*)");
    regex opt_upper_re("(.*)OPT(.*)");
    regex opt_ASEC_re("(.*)maxcyc=1,maxstep=15(.*)");
    regex espCharges("(.*)chelpg(.*)", regex_constants::icase);
    this->calcDone = 0;
    this->not_stop = 0;
    this->chelpg = 0;

    vector <string> fileLines;
    while(!arq.eof()){
        getline(arq, lineSTR);
        fileLines.push_back(lineSTR);
        if (regex_match(lineSTR, opt_re)){
            this->optAsw = 1;
        } else if (regex_match(lineSTR, opt_ASEC_re)){
            this->not_stop = 1;
        } else if (regex_match(lineSTR, opt_upper_re)){
            this->optAsw = 1;
        } else if (regex_match(lineSTR, states_re)){
            this->stateAsw = 1;
        } else if (regex_match(lineSTR, done_re)){
            this->calcDone = 1;
            vector <string> temp = splitString(lineSTR, ' ');
            this->date =  temp[6] + " " + temp[7] + " " + temp[8] + " " + temp[9] + " " + temp[10];
        } else if (regex_match(lineSTR, espCharges)){
            this->chelpg = 1;
        };
    };
    arq.close();
    if (!this->not_stop){
        this->molConstructor(fileLines);
        if (this->polarAsw){
            this->makePolar(fileLines);
        };
        if (this->stateAsw) {
            this->exSates = ExcStates(this->statesNum(fileLines));
            this->makeStates(fileLines);
        };
    } else{
        this->makeGradient(fileLines);
    };
    fileLines.clear();
    if (!this->not_stop){
        if (!this->calcDone){
            throw invalid_argument("The calculation, in this file " + filePath + ", has not been completed\n");
        };
    };
};

int G16LOGfile::statesNum(vector <string> fileLines){
    int statesNumber = 3;
    regex regexNstates1("(.*)TD(.*)nstates(.*)=(.*)[0-9]+(.*)");
    regex regexNstates2("(.*)TD(.*)NSTATES(.*)=(.*)[0-9]+(.*)");
    regex regexNstates3("(.*)TD(.*)Nstates(.*)=(.*)[0-9]+(.*)");
    for (int i = 0; i < (int) fileLines.size(); i++){
        if(regex_match(fileLines[i], regexNstates1) || regex_match(fileLines[i], regexNstates2) || regex_match(fileLines[i], regexNstates3)){
            if (statesNumber == 3){
                string temp = splitString(splitString(fileLines[i], '=').back(), ')')[0];
                statesNumber = stoi(temp);
            };
        };
    };
    return statesNumber;
};

void G16LOGfile::makeStates(vector <string> fileLines){
    regex excStatesRE(" Excited State(.*)[0-9]+:(.*)");
    bool takeLine = 0;
    int stateN = 0;
    vector <pair < pair <int, int>, double > > transitions;
    for (int i = 0; i < (int) fileLines.size(); i++){
        if(regex_match(fileLines[i], excStatesRE) && takeLine == 0){
            transitions.clear();
            takeLine = 1;
            vector <string> sLine = splitString(fileLines[i], ' ');
            string stateNtemp = sLine[2];
            stateNtemp.pop_back();
            stateN = stoi(stateNtemp);
            string sym = splitString(sLine[3], '-')[0];
            this->exSates.setSymmetry(stateN, sym);
            this->exSates.setEnergy(stateN, stod(sLine[4]));
            this->exSates.setWavelength(stateN, stod(sLine[6]));
            this->exSates.setOscillatorForce(stateN, stod(splitString(sLine[8], '=')[1]));
        } else if (takeLine == 1){
            vector <string> splitedLine = splitString(fileLines[i], ' ');
            if (splitedLine.size() != 3){
                this->exSates.setTransitions(stateN, transitions);
                takeLine = 0;
            } else {
                int trans1 = stoi(splitedLine[0]);
                int trans2 = stoi(splitString(splitedLine[1], '>')[1]);
                pair < pair <int, int>, double > trans = {pair <int, int> {trans1, trans2}, stod(splitedLine[2])};
                transitions.push_back(trans);
            };
        };
    };
};

void G16LOGfile::molConstructor(vector <string> fileLines){
    regex molecule_re1("(.*)Symbolic Z-matrix:");
    bool mol_re1_Bool = 0;
    regex molecule_re2("(.*)Standard orientation:(.*)");
    regex scf_re("(.*)SCF Done:(.*)");
    regex size_re("(.*)NAtoms=(.*)");
    regex omo_re("(.*)occ. eigenvalues(.*)");
    regex umo_re("(.*)virt. eigenvalues(.*)");
    regex basis_re(" Standard basis:(.*)");
    regex startCharge [1];
    int startMoleculeRef = 0;
    int chargeLine = 0;

    if (this->chelpg){
        regex chargeType("(.*)ESP charges:(.*)", regex_constants::icase);
        startCharge[0] = chargeType;
    } else {
        regex chargeType("(.*)Mulliken charges:(.*)", regex_constants::icase);
        startCharge[0] = chargeType;
    };
    for (int i = 0; i < (int) fileLines.size(); i++){
        if (regex_match(fileLines[i], scf_re)){
            vector <string> splittedLine = splitString(fileLines[i], ' ');
            this->energy = stod(splittedLine[4]);
            this->levelTheory = splitString(splitString(splittedLine[2], '(')[1], ')')[0];
        } else if (regex_match(fileLines[i], size_re)) {
            vector <string> splittedLine = splitString(fileLines[i], ' ');
            this->size = stoi(splittedLine[1]);
        } else if (regex_match(fileLines[i], omo_re)){
            vector <string> splittedLine = splitString(fileLines[i], ' ');
            for (int j = 4; j < (int) splittedLine.size(); j++){
                this->occOrb.push_back(stod(splittedLine[j]));
            };
        } else if (regex_match(fileLines[i], umo_re)){
            vector <string> splittedLine = splitString(fileLines[i], ' ');
            for (int j = 4; j < (int) splittedLine.size(); j++){
                this->virtOrb.push_back(stod(splittedLine[j]));
            };
        } else if (!this->optAsw){
            if (regex_match(fileLines[i], molecule_re1)) {
                startMoleculeRef = i+2;
                mol_re1_Bool = 1;
            };
        } else if (this->optAsw){
            if (regex_match(fileLines[i], molecule_re2)) {
                startMoleculeRef = i+5;
            };
        };
        if (this->basis.size() == 0){
            if (regex_match(fileLines[i], basis_re)){
                this->basis = splitString(fileLines[i], ' ')[2];
            };
        } else if (regex_match(fileLines[i], startCharge[0])){
            chargeLine = i + 2;
        };
        
    };
    int endCharge = chargeLine + this->size;
    vector <double> atomicCharges(this->size, 0.0);
    int count = 0;
    for (int i = chargeLine; i < endCharge; i++){
        vector <string> charString = splitString(fileLines[i], ' ');
        atomicCharges.at(count) = stod(charString.at(2));
        count ++;
    };

    int endMoleculeRef = startMoleculeRef + this->size;
    count = 0;
    for (int i = startMoleculeRef; i < endMoleculeRef; i++){
        vector <string> molLine = splitString(fileLines[i], ' ');
        if(mol_re1_Bool){
            this->molecule.addAtom(molLine[0], stod(molLine[1]), stod(molLine[2]), stod(molLine[3]), atomicCharges.at(count));
        } else {
            this->molecule.addAtom(molLine[1], stod(molLine[3]), stod(molLine[4]), stod(molLine[5]), atomicCharges.at(count));
        }
        count ++;
    };
};

void G16LOGfile::makePolar(vector <string> fileLines){
    regex dipole_re("(.*)Electric dipole moment(.*)input orientation(.*)");
    regex alpha_re("(.*)Dipole polarizability, Alpha(.*)input orientation(.*)");
    regex beta_re("(.*)First dipole hyperpolarizability, Beta(.*)input orientation(.*)");
    regex gamma_re("(.*)Second dipole hyperpolarizability, Gamma(.*)input orientation(.*)");
    int dipole_num_start = 0, dipole_num_end = 0;
    int alpha_num_start = 0, alpha_num_end = 0, beta_num_start = 0, beta_num_end = 0, gamma_num_start = 0, gamma_num_end = 0;
    for (int i = 0; i < (int) fileLines.size(); i++){
        if (regex_match(fileLines[i], dipole_re)){
            dipole_num_start = i+3;
            dipole_num_end = dipole_num_start + 4;
        } else if (regex_match(fileLines[i], alpha_re)){
            alpha_num_start = i+2;
        } else if (alpha_num_end == 0 && alpha_num_start != 0){
            if (fileLines[i] == ""){
                alpha_num_end = i;
            };
        } else if (regex_match(fileLines[i], beta_re)){
            beta_num_start = i+4;
        } else if (beta_num_end == 0 && beta_num_start != 0){
            if (fileLines[i] == ""){
                beta_num_end = i;
            };
        } else if (regex_match(fileLines[i], gamma_re)){
            gamma_num_start = i+4;
        } else if (gamma_num_end == 0 && gamma_num_start != 0){
            if (fileLines[i] == ""){
                gamma_num_end = i;
            };
        };
    };
    for (int i = dipole_num_start; i < dipole_num_end; i++){
        vector <string> dipLine = splitString(fileLines[i], ' ');
        vector <string> sValue = splitString(dipLine[2], 'D');
        this->polarValues.setDipole(dipLine[0], stod(sValue[0] + "e" + sValue[1]));
    };
    regex a_re("(.*)Alpha(.*)");
    string a = "";
    for (int i = alpha_num_start; i < alpha_num_end; i++){
        if (regex_match(fileLines[i], a_re)){
            a = fileLines[i];
        } else if(a != ""){
            vector <string> alpLine = splitString(fileLines[i], ' ');
            if (alpLine[2] != "esu)"){
                vector <string> aValue = splitString(alpLine[2], 'D');
                this->polarValues.setAlpha(splitString(a, ':')[0].erase(0, 1), alpLine[0], stod(aValue[0] + "e" + aValue[1]));
            };
        };
    };
    regex b_re("(.*)Beta(.*)");
    string b = "";
    for (int i = beta_num_start; i < beta_num_end; i++){
        if (regex_match(fileLines[i], b_re)){
            b = fileLines[i];
        } else if(b != ""){
            vector <string> betaLine = splitString(fileLines[i], ' ');
            if (betaLine[1] == "(z)"){
                vector <string> bValue = splitString(betaLine[3], 'D');
                this->polarValues.setBeta(splitString(b, ':')[0].erase(0, 1), "|| (z)", stod(bValue[0] + "e" + bValue[1]));
            } else if (betaLine[2] != "esu)"){
                vector <string> bValue = splitString(betaLine[2], 'D');
                this->polarValues.setBeta(splitString(b, ':')[0].erase(0, 1), betaLine[0], stod(bValue[0] + "e" + bValue[1]));
            };
        };
    };
    regex g_re("(.*)Gamma(.*)");
    string g = "";
    for (int i = gamma_num_start; i < gamma_num_end; i++){
        if (regex_match(fileLines[i], g_re)){
            g = fileLines[i];
        } else if(g != ""){
            vector <string> gammaLine = splitString(fileLines[i], ' ');
            if (gammaLine[2] != "esu)"){
                vector <string> gValue = splitString(gammaLine[2], 'D');
                this->polarValues.setGamma(splitString(g, ':')[0].erase(0, 1), gammaLine[0], stod(gValue[0] + "e" + gValue[1]));
            };
        };
    };
}

void G16LOGfile::makeGradient(vector <string> fileLines){
    regex grad_re("(.*)Forces \\(Hartrees/Bohr\\)(.*)");
    regex grad_re_end(" Cartesian Forces:  (.*)");
    int grad_num_start =0, grad_num_end = 0;
        for (int i = 0; i < (int) fileLines.size(); i++){
            if (regex_match(fileLines[i], grad_re)){
                grad_num_start = i+3;
            } else if (regex_match(fileLines[i], grad_re_end)){
                grad_num_end = i-1;
            };
        };
    vector <string> gradLine;
    vector < vector <double> > outsideVector;
    for (int i = grad_num_start; i < (int) grad_num_end; i++){
        gradLine.push_back(splitString(fileLines[i], ' ')[2]);
        gradLine.push_back(splitString(fileLines[i], ' ')[3]);
        gradLine.push_back(splitString(fileLines[i], ' ')[4]);
    };
    for (int i = 0; i < (int)gradLine.size(); i+=3){
        vector <double> tempGrad(3, 0);
        tempGrad.at(0) = stod(gradLine[i]);
        tempGrad.at(1) = stod(gradLine[i+1]);
        tempGrad.at(2) = stod(gradLine[i+2]);
        outsideVector.push_back(tempGrad); 
    };
    Matrix cGradient = Matrix(outsideVector);
    this->gradientValues.setGradient(cGradient);
};

Matrix G16LOGfile::getGradient(){
    return this->gradientValues.getGradient();

};
double G16LOGfile::scfEnergy(){
    return this->energy;
};

Molecule G16LOGfile::getMolecule(){
    return this->molecule;
};

double G16LOGfile::getDipole(string name){
    if (this->polarAsw){
        return this->polarValues.getDipole(name);
    } else {
        return 0.0;
    };
};

double G16LOGfile::getAlpha(string eleName, string name){
    if (this->polarAsw){
        return this->polarValues.getAlpha(eleName, name);
    } else {
        return 0.0;
    };
};

double G16LOGfile::getBeta(string eleName, string name){
    if (this->polarAsw){
        return this->polarValues.getBeta(eleName, name);
    } else {
        return 0.0;
    };
};

double G16LOGfile::getGamma(string eleName, string name){
    if (this->polarAsw){
        return this->polarValues.getGamma(eleName, name);
    } else {
        return 0.0;
    };
};

double G16LOGfile::getOscillatorForce(int state){
    double result = 0.0;
    if (this->stateAsw){
        result =  this->exSates.getOscillatorForce(state);
    };
    return result;
};

double G16LOGfile::getWavelength(int state){
    double result = 0.0;
    if (this->stateAsw){
        result = this->exSates.getWavelength(state);
    };
    return result;
};

string G16LOGfile::getSymmetry(int state){
    string result = "";
    if (this->stateAsw){
        result = this->exSates.getSymmetry(state);
    };
    return result;
};

vector <double> G16LOGfile::getOscillatorForces(){
    vector <double> result (1, 0.0);
    if (this->stateAsw){
        result.clear();
        for (int i = 1; i < this->exSates.getstatesNumber()+1; i++){
            result.push_back(this->exSates.getOscillatorForce(i));
        };
    };
    return result;
};

vector <double> G16LOGfile::getWavelengths(){
    vector <double> result (1, 0.0);
    if (this->stateAsw){
        result.clear();
        for (int i = 1; i < this->exSates.getstatesNumber()+1; i++){
            result.push_back(this->exSates.getWavelength(i));
        };
    };
    return result;
};

vector <string> G16LOGfile::getSymmetries(){
    vector <string> result (1, "");
    if (this->stateAsw){
        result.clear();
        for (int i = 1; i < this->exSates.getstatesNumber()+1; i++){
            result.push_back(this->exSates.getSymmetry(i));
        };
    };
    return result;
};

vector <string> G16LOGfile::getTransition(int state){
    vector <string> result(1, "");
    if (this->stateAsw){
        result.clear();
        vector <pair < pair <int, int>, double > > temp = this->exSates.getTransition(state);
        for (int i = 0; i< (int) temp.size(); i++){
            string temp1 = to_string(temp[i].first.first) + " -> " + to_string(temp[i].first.second) + "   " + to_string(temp[i].second);
            result.push_back(temp1);
        };
    };
    return result;
};

vector < vector <string>> G16LOGfile::getTransitions(){
    vector < vector <string>> result (1);
    if (this->stateAsw){
        result.clear();
        for (int i = 1; i < this->exSates.getstatesNumber()+1; i++){
            vector <string> temp = this->getTransition(i);
            result.push_back(temp);
        };
    };
    return result;
};

vector <string> G16LOGfile::getTransitionsStr(){
    vector <string> results (1, "");
    if (this->stateAsw){
        results.clear();
        string result = "";
        for (int j = 1; j < this->exSates.getstatesNumber() + 1; j++){
            vector <pair < pair <int, int>, double > > temp = this->exSates.getTransition(j);
            for (int i = 0; i < (int) temp.size(); i++){
                int firstOrb = temp[i].first.first;
                int secondOrb = temp[i].first.second;
                double coef = temp[i].second;
                string t =  "Orb. " + to_string(firstOrb) + " to Orb. " + to_string(secondOrb) + ", L.C. coef.: " + to_string(coef);
                if (result.size() == 0){
                    result = t;
                } else {
                    result = result + "\n" + t;
                };
            };
            results.push_back(result);
        }
    };
    return results;
};

string G16LOGfile::toStr(){
    string result = this->molecule.toStr();
    result = result + "\nEnergy: " + to_string(this->energy) + " hartree";
    result = result + "\nCalculation done in " + this->date + " Level of Theory: " + this->levelTheory + "/" + this->basis;
    return result;
};

vector <vector <string>> G16LOGfile::getTransContributions(){
    vector <vector <string>> results (1);
    if (this->stateAsw){
        results.clear();
        vector <string> result;
        for (int j = 1; j < this->exSates.getstatesNumber() + 1; j++){
            vector < pair <string, double> > temp = this->exSates.getTransContribution(j);
            for (int i = 0; i < (int) temp.size(); i++){
                string orbs = temp[i].first;
                double cont = temp[i].second;
                string t =  orbs + ", Contrib.: " + to_string(cont);
                result.push_back(t);
            };
            results.push_back(result);
        }
    };
    return results;
};





/*
----------------------- G16FCHKfile -----------------------
*/

G16FCHKfile::G16FCHKfile(string filePath){
    ifstream arq;
    arq.open(filePath, ifstream::in);
    string lineSTR;


    vector <string> fileLines;
    while(!arq.eof()){
        getline(arq, lineSTR);
        fileLines.push_back(lineSTR);
    };
    arq.close();

    this->molConstructor(fileLines);
    this->makeGradient(fileLines);
};

void G16FCHKfile::molConstructor(vector <string> fileLines){
    
    regex molecule_atoms_re("(.*)Atomic numbers(.*)");
    regex molecule_atoms_re_end("(.*)Nuclear charges(.*)");
    regex molecule_re("(.*)Current cartesian coordinates(.*)");
    regex molecule_re_end("(.*)Number of symbols in(.*)");
    regex charge_type("(.*)chelpg(.*)", regex_constants::icase);
    regex save_charges("(.*)SaveESP(.*)", regex_constants::icase);
    regex molecule_charge [1];
    regex scf_re("(.*)SCF Energy(.*)");
    regex size_re("(.*)Number of atoms(.*)");
    regex ele_re("(.*)Number of electrons(.*)");
    int startMoleculeRef = 0;
    int endMoleculeRef = 0;
    int startAtomsRef = 0;
    int endAtomsRef = 0;
    int startCharges = 0;
    int endCharge = 0;
    regex chargeStringMul("(.*)Mulliken Charges(.*)", regex_constants::icase);
    molecule_charge[0] = chargeStringMul;

    
    for (int i = 0; i < (int) fileLines.size(); i++){
        if (regex_match(fileLines[i], charge_type)) {
            if (regex_match(fileLines[i], save_charges)) {
                regex chargeString("(.*)MM charges(.*)", regex_constants::icase);
                molecule_charge[0] = chargeString;
            } else {
                regex chargeString("(.*)ESP Charges(.*)", regex_constants::icase);
                molecule_charge[0] = chargeString;

            };
        };
    };
    for (int i = 0; i < (int) fileLines.size(); i++){
        if (regex_match(fileLines[i], scf_re)){
            vector <string> splittedLine = splitString(fileLines[i], ' ');
            this->energy = stod(splittedLine[3]);
        } else if (regex_match(fileLines[i], ele_re)) {
            vector <string> splittedLine = splitString(fileLines[i], ' ');
            this->electronNumber = stoi(splittedLine[4]);
        } else if (regex_match(fileLines[i], size_re)) {
            vector <string> splittedLine = splitString(fileLines[i], ' ');
            this->size = stoi(splittedLine[4]);
        } else if (regex_match(fileLines[i], molecule_atoms_re)) {
                startAtomsRef = i+1;
        } else if (regex_match(fileLines[i], molecule_atoms_re_end)) {
                endAtomsRef = i;
        } else if (regex_match(fileLines[i], molecule_re)) {
                startMoleculeRef = i + 1;
        } else if (regex_match(fileLines[i], molecule_re_end)) {
                endMoleculeRef = i;
        } else if (regex_match(fileLines[i], molecule_charge[0])) {
                startCharges = i + 1;
        };
        
    };
    vector <double> atomsVector;
    for (int i = startAtomsRef; i < endAtomsRef; i++){
        for (int j = 0; j < (int) splitString(fileLines[i], ' ').size(); j++){
            atomsVector.push_back(stoi(splitString(fileLines[i], ' ')[j]));
        };
    };
    if (atomsVector.size() % 5 == 0){
        endCharge = startCharges + atomsVector.size()/5; 
    } else {
                endCharge = startCharges + (atomsVector.size()/5)+1; 

    };
    vector <double> atomsCharge(atomsVector.size(), 0.0);
    int c = 0;
    for (int i = startCharges; i < endCharge; i++){
        for (int j = 0; j < (int) splitString(fileLines[i], ' ').size(); j++){
            atomsCharge.at(c) = (stod(splitString(fileLines[i], ' ')[j]));
            c++;
        };
    };
    vector <double> molLine; 
    vector < vector <double> > molVector;
    for (int i = startMoleculeRef; i < (int) endMoleculeRef; i++){
        for (int j = 0; j < (int) splitString(fileLines[i], ' ').size(); j++){
            molLine.push_back(stod(splitString(fileLines[i], ' ')[j]) * 0.529177249);
        };
    };
    for (int i = 0; i < (int) molLine.size(); i+=3){
        vector <double> tempGrad(3, 0); 
        tempGrad.at(0) = molLine[i];
        tempGrad.at(1) = molLine[i+1];
        tempGrad.at(2) = molLine[i+2];
        molVector.push_back(tempGrad); 
    };
    for (int i = 0; i < (int) atomsVector.size(); i++){
        this->molecule.addAtom(atomsVector[i], molVector[i][0], molVector[i][1], molVector[i][2], atomsCharge[i]);
    };
    
};

void G16FCHKfile::makeGradient(vector <string> fileLines){
    regex start_gradient_re("Cartesian Gradient (.*)");
    regex end_gradient_re("Nonadiabatic coupling (.*)");
    int start = 0, end = 0;
    for (int i = 0; i < (int) fileLines.size(); i++){
        if (regex_match(fileLines[i], start_gradient_re)){
            start = i+1;
        } else if (regex_match(fileLines[i], end_gradient_re)){
            end = i;
        };
    };
    vector <string> gradLine; 
    vector < vector <double> > outsideVector;
    for (int i = start; i < (int) end; i++){
        for (int j = 0; j < (int) splitString(fileLines[i], ' ').size(); j++){
            gradLine.push_back(splitString(fileLines[i], ' ')[j]);
        };
    };
    for (int i = 0; i < (int) gradLine.size(); i+=3){
        vector <double> tempGrad(3, 0); 
        tempGrad.at(0) = stod(gradLine[i]);
        tempGrad.at(1) = stod(gradLine[i+1]);
        tempGrad.at(2) = stod(gradLine[i+2]);
        outsideVector.push_back(tempGrad); 
    };
    Matrix cGradient = Matrix(outsideVector);
    this->gradientValues.setGradient(cGradient);
};

Molecule G16FCHKfile::getMolecule(){
    return this->molecule;
};

Matrix G16FCHKfile::getCartesianGradient(){
    return this->gradientValues.getGradient();
};

/*
----------------------- GradientValues -----------------------
*/

GradientValues::GradientValues(){
};

void GradientValues::setGradient(Matrix gradient){
    this->cartesianGradient = gradient;
};

Matrix GradientValues::getGradient(){
    return this->cartesianGradient;
};

/*
----------------------- PolarValues -----------------------
*/

PolarValues::PolarValues(){
};

void PolarValues::setDipole(string name, double value){
    this->dName.push_back(name);
    this->dValue.push_back(value);
};

void PolarValues::setAlpha(string eleName, string name, double value){
    int n = -1;
    for (int i = 0; i < (int) this->aName.size(); i++){
        if (eleName == aName[i].first){
            n = i;
        };
    };
    if (n != -1){
        vector<string> tN = this->aName.at(n).second;
        vector<double> tV = this->aValue.at(n).second;
        tN.push_back(name);
        tV.push_back(value);
        pair< string, vector<string>> tempN = {eleName, tN};
        pair< string, vector<double>> tempV = {eleName, tV};
        this->aName.at(n) = tempN;
        this->aValue.at(n) = tempV;
    } else if (n == -1){
        pair< string, vector<string>> tempN = {eleName, vector<string>(1, name)};
        pair< string, vector<double>> tempV = {eleName, vector<double>(1, value)};
        this->aName.push_back(tempN);
        this->aValue.push_back(tempV);
    };
};

void PolarValues::setBeta(string eleName, string name, double value){
    int n = -1;
    for (int i = 0; i < (int) this->bName.size(); i++){
        if (eleName == bName[i].first){
            n = i;
        };
    };
    if (n != -1){
        vector<string> tN = this->bName.at(n).second;
        vector<double> tV = this->bValue.at(n).second;
        tN.push_back(name);
        tV.push_back(value);
        pair< string, vector<string>> tempN = {eleName, tN};
        pair< string, vector<double>> tempV = {eleName, tV};
        this->bName.at(n) = tempN;
        this->bValue.at(n) = tempV;
    } else if (n == -1){
        pair< string, vector<string>> tempN = {eleName, vector<string>(1, name)};
        pair< string, vector<double>> tempV = {eleName, vector<double>(1, value)};
        this->bName.push_back(tempN);
        this->bValue.push_back(tempV);
    };
};

void PolarValues::setGamma(string eleName, string name, double value){
    int n = -1;
    for (int i = 0; i < (int) this->gName.size(); i++){
        if (eleName == gName[i].first){
            n = i;
        };
    };
    if (n != -1){
        vector<string> tN = this->gName.at(n).second;
        vector<double> tV = this->gValue.at(n).second;
        tN.push_back(name);
        tV.push_back(value);
        pair< string, vector<string>> tempN = {eleName, tN};
        pair< string, vector<double>> tempV = {eleName, tV};
        this->gName.at(n) = tempN;
        this->gValue.at(n) = tempV;
    } else if (n == -1){
        pair< string, vector<string>> tempN = {eleName, vector<string>(1, name)};
        pair< string, vector<double>> tempV = {eleName, vector<double>(1, value)};
        this->gName.push_back(tempN);
        this->gValue.push_back(tempV);
    };
};

double PolarValues::getDipole(string name){
    int n = 0;
    for (int i = 0; i < (int) this->dName.size(); i++){
        if (this->dName[i] == name){
            n = i;
        };
    };
    return this->dValue[n];
};

double PolarValues::getAlpha(string eleName, string name){
    int n = 0, m = 0;
    for (int i = 0; i < (int) this->aName.size(); i++){
        if (this->aName[i].first == eleName){
            n = i;
        };
    };
    for (int i = 0; i < (int) this->aName[i].second.size(); i++){
        if (this->aName[n].second[i] == name){
            m = i;
        };
    };
    return this->aValue[n].second[m];
};

double PolarValues::getBeta(string eleName, string name){
     int n = 0, m = 0;
    for (int i = 0; i < (int) this->bName.size(); i++){
        if (this->bName[i].first == eleName){
            n = i;
        };
    };
    for (int i = 0; i < (int) this->bName[i].second.size(); i++){
        if (this->bName[n].second[i] == name){
            m = i;
        };
    }
    return this->bValue[n].second[m];
};

double PolarValues::getGamma(string eleName, string name){
     int n = 0, m = 0;
     for (int i = 0; i < (int) this->gName.size(); i++){
         if (this->gName[i].first == eleName){
             n = i;
         };
     };
     for (int i = 0; i < (int) this->gName[i].second.size(); i++){
         if (this->gName[n].second[i] == name){
             m = i;
         };
     }
     return this->gValue[n].second[m];
};

/*
----------------------- ExcStates -----------------------
*/
ExcStates::ExcStates(){
    this->statesNumber = 0;
};

ExcStates::ExcStates(int statesNumber){
    this->statesNumber = statesNumber;
    this->wlValues.resize(this->statesNumber);
    this->energies.resize(this->statesNumber);
    this->oscillator.resize(this->statesNumber);
    this->symmetries.resize(this->statesNumber);
    this->transitions.resize(this->statesNumber);
};

void ExcStates::setWavelength(int state, double value){
    this->wlValues.at(state-1) = value;
};

void ExcStates::setEnergy(int state, double value){
    this->energies.at(state-1) = value;
};

void ExcStates::setOscillatorForce(int state, double value){
    this->oscillator.at(state-1) = value;
};

void ExcStates::setSymmetry(int state, string value){
    this->symmetries.at(state-1) = value;
};

void ExcStates::setTransitions(int state, vector <pair < pair <int, int>, double > > values){
    this->transitions.at(state-1) = values;
};

string ExcStates::getSymmetry(int state){
    return this->symmetries[state-1];
};

double ExcStates::getWavelength(int state){
    return this->wlValues[state-1];
};
double ExcStates::getEnergy(int state){
    return this->energies[state-1];
};
double ExcStates::getOscillatorForce(int state){
    return this->oscillator[state-1];
};

vector <pair < pair <int, int>, double > > ExcStates::getTransition(int state){
    return this->transitions[state-1];
};

vector < pair <string, double> > ExcStates::getTransContribution(int state){
    vector <pair < pair <int, int>, double > > trans = this->transitions[state-1];
    vector < pair <string, double> > result;
    for (int i = 0; i < (int) trans.size(); i++){
        string orbTrans;
        orbTrans = to_string(trans[i].first.first) + " -> " + to_string(trans[i].first.second);
        double contrib = 2 * pow(trans[i].second, 2);
        result.push_back(pair <string, double> {orbTrans, contrib});
    };
    return result;
};

int ExcStates::getstatesNumber(){
    return this->statesNumber;
};
