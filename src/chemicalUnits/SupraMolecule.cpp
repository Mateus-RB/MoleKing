//   MoleKing //
//
//   File:        [SupraMolecule.cpp]
//
//   Author(s):   ['LEEDMOL Reserch Group']
//   Site(s):     ['https://www.researchgate.net/lab/LEEDMOL-Heibbe-Cristhians-Research-Group-Heibbe-Cristhian-B-De-Oliveira']
//   Email(s):    ['heibbe@ufg.br']
//   Credits:     ['Copyright © 2023 LEEDMOL. All rights reserved.']
//   Date:        ['17.01.2023']
//   Version:     ['1.5.1']
//   Status:      ['Development']
//   Language:    ['C++','Python']
//   Description: ['A python module written in C++ for theoretical chemestry']


#include "SupraMolecule.hpp"

typedef vector <vector <int> > VectorsInt;

SupraMolecule::SupraMolecule(){};

void SupraMolecule::setCharge(){
    this->charge = 0;
    for (int i = 0; i < (int) this->supraMolecule.size(); i++){
        this->charge = this->charge + supraMolecule[i].getCharge();
    };
};

void SupraMolecule::addMolecule(Molecule molecule){
    this->supraMolecule.push_back(molecule);
    this->setCharge();
};

void SupraMolecule::addAtomToMolecule(int molNumber, Atom atom){
    this->supraMolecule[molNumber].addAtom(atom);
    this->setCharge();
};

void SupraMolecule::addAtomToMolecule(int molNumber, string atomSymbol, double xPos, double yPos, double zPos){
    this->supraMolecule[molNumber].addAtom(atomSymbol, xPos, yPos, zPos);
    this->setCharge();
};

string SupraMolecule::toStr(){
    string result = "Supramolecule: ";

    // result = result << this->supraMolecule.size() << " molecules with a total charge of " << this->charge << endl;

    result = result + to_string(this->supraMolecule.size()) + " molecules with a total charge of " + to_string(this->charge) + ".";

    for (int i = 0; i < (int) this->supraMolecule.size(); i++){
        result = result + "\n    " + this->supraMolecule[i].toStr();
    };
    return result;
}

Molecule SupraMolecule::getMolecule(int numberMolecule){
    return this->supraMolecule[numberMolecule];
};

void SupraMolecule::setMultiplicity(int multiplicity){
    this->multiplicity = multiplicity;
};

int SupraMolecule::getMultiplicity(){
    return this->multiplicity;
};

double SupraMolecule::getCharge(){
    return this->charge;
};

long SupraMolecule::getSize(){
    return this->supraMolecule.size();
};

Point SupraMolecule::getMassCenter(){
    vector <double> massVector;
    vector <double> coordX;
    vector <double> coordY;
    vector <double> coordZ;
    for (int i = 0; i < (int) this->supraMolecule.size(); i++){
        for (int j = 0; j < this->supraMolecule[i].getSize(); j++){
            Atom atom = this->supraMolecule[i].getAtomObj(j);
            massVector.push_back(atom.getAtomicMass());
            coordX.push_back(atom.getX());
            coordY.push_back(atom.getY());
            coordZ.push_back(atom.getZ());
        };
    };
    return MassCenter(massVector, coordX, coordY, coordZ).getMassCenter();
};

void SupraMolecule::translation(Vector3D traslationVector){
    for(int i = 0; i < (int) this->supraMolecule.size(); i++){
        this->supraMolecule[i].translation(traslationVector);
    };
};

void SupraMolecule::moveMassCenter(double x, double y, double z){
    Vector3D traslationVector = Vector3D({x, y, z}, this->getMassCenter().getCoords('c'));
    this->translation(traslationVector);
};

void SupraMolecule::moveTail(int molNumber, int atomNumber, double x, double y, double z){
    vector <double> pos = this->supraMolecule.at(molNumber).getAtomObj(atomNumber).getPos();
    Vector3D traslationVector = Vector3D({x, y, z}, pos);
    this->translation(traslationVector);
};

void SupraMolecule::spinSupraMolecule(double angle, char axis){
    if (axis == 'x') {
        Vector3D spinVector = Vector3D({1.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
        this->spinSupraMolecule(angle , spinVector);
    } else if (axis == 'y'){
        Vector3D spinVector = Vector3D({0.0, 1.0, 0.0}, {0.0, 0.0, 0.0});
        this->spinSupraMolecule(angle , spinVector);
    } else {
        Vector3D spinVector = Vector3D({0.0, 0.0, 1.0}, {0.0, 0.0, 0.0});
        this->spinSupraMolecule(angle , spinVector);
    };
};

void SupraMolecule::spinSupraMolecule(double angle, Vector3D spinVector){
    for (int i = 0; i < (int) this->supraMolecule.size(); i++){
        for (int j = 0; j < (int) this->supraMolecule[i].getSize(); j++){
            this->supraMolecule[i].getAtomObj(j).rotationAxis(angle, spinVector);
        };
    };
};


void SupraMolecule::getMoleculeBonds(){
    bonds.clear();
    for (int i = 0; i < (int) this->supraMolecule.size(); i++){
        VectorsInt temp = this->supraMolecule[i].getIRCBonds();
        bonds.push_back(temp);
    };
};

void SupraMolecule::getMoleculeAngles(){
    angles.clear();
    for (int i = 0; i < (int) this->supraMolecule.size(); i++){
        VectorsInt temp = this->supraMolecule[i].getIRCAngles();
        bonds.push_back(temp);
    };
};

void SupraMolecule::getMoleculeTorsions(){
    dihedrals.clear();
    for (int i = 0; i < (int) this->supraMolecule.size(); i++){
        VectorsInt temp = this->supraMolecule[i].getIRCDihedrals();
        bonds.push_back(temp);
    };
};

// void SupraMolecule::standardOrientation(int molNumber){
//     vector <int> molAxis = this->supraMolecule[molNumber].molecularAxis();
//     Vector3D moveTailVec = Vector3D({0, 0, 0}, this->supraMolecule[molNumber].getAtomObj(molAxis[0]).getPos());
//     for (int i = 0; i < (int) this->supraMolecule.size(); i++){
//         this->supraMolecule[i].translation(moveTailVec);
//     };
//     vector <double> angles = this->supraMolecule[molNumber].standardOrientationPath();
//     Vector3D mMassCent = Vector3D({0, 0, 0}, this->supraMolecule[molNumber].getMassCenter().getCoords('c'));
//     this->supraMolecule[molNumber].translation(mMassCent);
//     for (int i = 0; i < (int) this->supraMolecule.size(); i++){
//         if (i != molNumber){
//             this->supraMolecule[i].spinMolecule(angles[0], 'y');
//             this->supraMolecule[i].spinMolecule(angles[1], 'z');
//             this->supraMolecule[i].spinMolecule(angles[2], 'x');
//             this->supraMolecule[i].spinMolecule(90, 'z');
//             this->supraMolecule[i].spinMolecule(90, 'y');
//             this->supraMolecule[i].translation(mMassCent);
//         };
//     };
// };

vector <VectorsInt> SupraMolecule::getIRCBonds(){
    vector <VectorsInt> irc;
    for (int i = 0; i < (int) this->supraMolecule.size(); i++){
        irc.push_back(this->supraMolecule[i].getIRCBonds());
    }
    return irc;
};

vector <VectorsInt> SupraMolecule::getIRCAngles(){
    vector <VectorsInt> irc;
    for (int i = 0; i < (int) this->supraMolecule.size(); i++){
        irc.push_back(this->supraMolecule[i].getIRCAngles());
    }
    return irc;
};

vector <VectorsInt> SupraMolecule::getIRCDihedrals(){
    vector <VectorsInt> irc;
    for (int i = 0; i < (int) this->supraMolecule.size(); i++){
        irc.push_back(this->supraMolecule[i].getIRCDihedrals());
    }
    return irc;
};

void SupraMolecule::removeAtom(int molNumber, int atomNumber){
    this->supraMolecule[molNumber].removeAtom(atomNumber);
};

void SupraMolecule::removeAtom(Atom atom){
    for (int i = 0; i < (int) this->supraMolecule.size(); i++){
        this->supraMolecule[i].removeAtom(atom);
    };
};

void SupraMolecule::removeElement(int molNumber, string element){
    this->supraMolecule[molNumber].removeElement(element);
};

void SupraMolecule::removeElement(string element){
    for (int i = 0; i < (int) this->supraMolecule.size(); i++){
        this->supraMolecule[i].removeElement(element);
    };
};


void SupraMolecule::removeMolecule(int molNumber){
    this->supraMolecule.erase(this->supraMolecule.begin()+molNumber);
};

void SupraMolecule::removeMolecule(Molecule molecule){
    for (int i = 0; i < (int) this->supraMolecule.size(); i++){
        if (this->supraMolecule[i] == molecule){
            this->supraMolecule.erase(this->supraMolecule.begin()+i);
        };
    };
};

double SupraMolecule::getSupraMolecularMass(){
    double supraMolecularMass = 0;
    for (int i = 0; i < (int) this->supraMolecule.size(); i++){
        supraMolecularMass = supraMolecularMass + supraMolecule[i].getMolecularMass();
    };
    return supraMolecularMass;
};

vector <Molecule> SupraMolecule::moleculeList(){
    return this->supraMolecule;
};

SupraMolecule::iterator SupraMolecule::begin(){
    return this->supraMolecule.begin();
};

SupraMolecule::iterator SupraMolecule::end(){
    return this->supraMolecule.end();
};

bool SupraMolecule::operator==(SupraMolecule sMol){
    if ((int) this->supraMolecule.size() == (int) sMol.getSize()){
        for (int i = 0; i < (int) this->supraMolecule.size(); i++){
            if (this->supraMolecule[i] == sMol.getMolecule(i)){
                continue;
            } else {
                return 0;
            };
        };
    } else{
        return 0;
    };
    return 1;
};

bool SupraMolecule::operator!=(SupraMolecule sMol){
    if ((int) this->supraMolecule.size() == (int) sMol.getSize()){
        for (int i = 0; i < (int) this->supraMolecule.size(); i++){
            if (this->supraMolecule[i] == sMol.getMolecule(i)){
                continue;
            } else {
                return 1;
            };
        };
    } else{
        return 1;
    };
    return 0;
};

bool SupraMolecule::operator<(SupraMolecule sMol){
    if ((int) this->supraMolecule.size() < (int) sMol.getSize()){
        return 1;
    } else {
        for (int i = 0; i < (int) this->supraMolecule.size(); i++){
            if (this->supraMolecule[i] < sMol.getMolecule(i)){
                continue;
            } else {
                return 0;
            };
        };
    };
};

bool SupraMolecule::operator>(SupraMolecule sMol){
    if ((int) this->supraMolecule.size() > (int) sMol.getSize()){
        return 1;
    } else {
        for (int i = 0; i < (int) this->supraMolecule.size(); i++){
            if (this->supraMolecule[i] > sMol.getMolecule(i)){
                continue;
            } else {
                return 0;
            };
        };
    };
};