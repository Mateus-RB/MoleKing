//   MoleKing //
//
//   File:        [pov.cpp]
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

#include "pov.hpp"


PovRay::PovRay(const Molecule& mol) : mol(mol) {
    this->mol = mol;
    this->PT = PeriodicTable(); 
    getLookAtPos();
    atomFarFromCM();
    buildFile();


    //todo: bond tickness, bond color, spin molecule, set povray name file, change atom radii, add vdw and cpk (to vdw just multiply by vdw radii [1.3])
    

    //! transparent render: povray molecule.pov +ua

}

string PovRay::getCamPos(){
    return "<" + to_string(this->mcCoords[0]+this->maxDist+1) + ", " + to_string(this->mcCoords[1]+this->maxDist+1) + ", " + to_string(this->mcCoords[2]+this->maxDist+1) + ">";
}

string PovRay::getLookAtPos(){
    this->mc = this->mol.getMassCenter();
    this->mcCoords = mc.getCoords('c');
    return "<" + to_string(this->mcCoords[0]) + ", " + to_string(this->mcCoords[1]) + ", " + to_string(this->mcCoords[2]) + ">";
}

void PovRay::atomFarFromCM(){

    double dist = 0;
    double maxDist = 0;
    double maxX = 0;
    double maxY = 0;
    double maxZ = 0;

    for (int i = 0; i < this->mol.getSize(); i++)
    {
        double cmX = this->mcCoords[0];
        double cmY = this->mcCoords[1];
        double cmZ = this->mcCoords[2];

        double atomX = this->mol.getAtomObj(i).getX();
        double atomY = this->mol.getAtomObj(i).getY();
        double atomZ = this->mol.getAtomObj(i).getZ();

        dist = sqrt(pow(cmX - atomX, 2) + pow(cmY - atomY, 2) + pow(cmZ - atomZ, 2));

        if (dist > maxDist)
        {
            maxDist = dist;
        }
    }
    this->maxDist = maxDist;
}

void PovRay::buildFile(){

    ofstream povFile;
    povFile.open("molecule.pov");
    povFile << "#declare BSAMBI = 0.2;\n#declare BSDIFF = 0.8;\n#declare BSSPEC = 0.3;\n#declare TRANS=0.0;\n\n"
            << "camera {\n"
            << "perspective\n"
            //<< "location <0, 0, 10>\n"
            << "location " << getCamPos() << "\n"
            << "look_at " << getLookAtPos() << "\n"
            << "}\n"
            << "background { rgbt " << HEX2RGB("FFFFFF") << ", 1> }\n"
            << "light_source { <500.0, 500.0, 500.0> rgb " << HEX2RGB("FFFFFF") << "> }\n"
            << "plane {<1.0, 1.0, 1.0>, 0.0\n"
            <<  "pigment{ color rgbt <1,1,1,1> }\n"
            << "finish {ambient BSAMBI diffuse BSDIFF specular BSSPEC }\n"
            << "no_shadow\n"
            << "}\n"
            << "#declare molecule = union {\n";
            for (int i = 0; i < this->mol.getSize(); i++)
            {
                povFile << "    sphere { <" << setAtomicPosition(i)[0] << ", " << setAtomicPosition(i)[1] << ", " << setAtomicPosition(i)[2] << ">  " << setAtomicRadii(this->mol.getAtomObj(i).getAtomicSymbol()) << "\n"
                << "        no_shadow\n"
                << "        pigment { color rgbt" << HEX2RGB(setAtomicColor(this->mol.getAtomObj(i).getAtomicSymbol())) << ", TRANS> }\n"
                << "        finish {ambient BSAMBI diffuse BSDIFF specular BSSPEC}\n"
                << "    }\n";
            }
    povFile << "}\n"
            << "object { molecule }\n";
    povFile.close();
}

vector <double> PovRay::setAtomicPosition(int atomicNumber){
    return {this->mol.getAtomObj(atomicNumber).getX(), this->mol.getAtomObj(atomicNumber).getY(), this->mol.getAtomObj(atomicNumber).getZ()};
}

double PovRay::setAtomicRadii(string atomicSymbol){
    return this->PT.getCovalentRadii(atomicSymbol);
}

string PovRay::setAtomicColor(string atomicSymbol){
    return this->PT.getColor(atomicSymbol);
}

string PovRay::HEX2RGB(string hex)
{
    map<char, int> hex2int;
    hex2int.insert(pair<char, int>('0', 0));
    hex2int.insert(pair<char, int>('1', 1));
    hex2int.insert(pair<char, int>('2', 2));
    hex2int.insert(pair<char, int>('3', 3));
    hex2int.insert(pair<char, int>('4', 4));
    hex2int.insert(pair<char, int>('5', 5));
    hex2int.insert(pair<char, int>('6', 6));
    hex2int.insert(pair<char, int>('7', 7));
    hex2int.insert(pair<char, int>('8', 8));
    hex2int.insert(pair<char, int>('9', 9));
    hex2int.insert(pair<char, int>('A', 10));
    hex2int.insert(pair<char, int>('B', 11));
    hex2int.insert(pair<char, int>('C', 12));
    hex2int.insert(pair<char, int>('D', 13));
    hex2int.insert(pair<char, int>('E', 14));
    hex2int.insert(pair<char, int>('F', 15));
    hex2int.insert(pair<char, int>('a', 10));
    hex2int.insert(pair<char, int>('b', 11));
    hex2int.insert(pair<char, int>('c', 12));
    hex2int.insert(pair<char, int>('d', 13));
    hex2int.insert(pair<char, int>('e', 14));
    hex2int.insert(pair<char, int>('f', 15));

    double r = hex2int[hex[0]] * 16 + hex2int[hex[1]];
    double g = hex2int[hex[2]] * 16 + hex2int[hex[3]];
    double b = hex2int[hex[4]] * 16 + hex2int[hex[5]];

    r = r / 255.0;
    g = g / 255.0;
    b = b / 255.0;
     
    return "<" + to_string(r) + ", " + to_string(g) + ", " + to_string(b);
}