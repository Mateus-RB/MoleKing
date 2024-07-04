//   MoleKing //
//
//   File:        [Molecule.cpp]
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

#include "Molecule.hpp"


double Molecule::angleToSpinInAref(int ref, char axisName){
    vector <double> cart = this->molecule[ref].getPos();
    if (axisName == 'x'){
        Point victor = Point(cart[0], cart[1], cart[2], 'c');
        Point victorDoidera = Point(cart[0], cart[1], cart[2], 'c');
        Vector3D xAxis = Vector3D({1.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
        victorDoidera.rotationVector(180, xAxis);
        Vector3D victores = Vector3D(victor.getCoords('c'), victorDoidera.getCoords('c'));
        double linhaDaLoucura = victores.magnitude();
        double raioVictoral = linhaDaLoucura/2;
        double xDaLoucura = -1 *(sqrt(pow(raioVictoral, 2) - pow(victorDoidera.getCoords('c')[2],2)) - victorDoidera.getCoords('c')[0]);
        Point centroDaLoucura = Point(xDaLoucura, victor.getCoords('c')[2], 0.0, 'c');
        Vector3D  sonOfVictor = Vector3D(victor.getCoords('c'), centroDaLoucura.getCoords('c'));
        double xTombado = raioVictoral + xDaLoucura;
        Vector3D sonOfZ = Vector3D( {xTombado , victor.getCoords('c')[1], 0.0} , centroDaLoucura.getCoords('c'));
        double anguloDaLoucura = sonOfVictor.angle(sonOfZ);
        return anguloDaLoucura;
    } 
    else 
    {
        Point victor = Point(cart[0], cart[1], cart[2], 'c');
        Point victorDoidera = Point(cart[0], cart[1], cart[2], 'c');
        Vector3D yAxis = Vector3D({0.0, 1.0, 0.0}, {0.0, 0.0, 0.0});
        victorDoidera.rotationVector(180, yAxis);
        // cout << "Victor: " << victor.getCoords('c')[0] << " " << victor.getCoords('c')[1] << " " << victor.getCoords('c')[2] << endl;
        // cout << "VictorDoidera: " << victorDoidera.getCoords('c')[0] << " " << victorDoidera.getCoords('c')[1] << " " << victorDoidera.getCoords('c')[2] << endl;

        Vector3D victores = Vector3D(victor.getCoords('c'), victorDoidera.getCoords('c'));
        double linhaDaLoucura = victores.magnitude();
        // cout << "Linha da Loucura: " << linhaDaLoucura << endl;

        double raioVictoral = linhaDaLoucura/2;
        cout << "RaioVictoral: " << raioVictoral << endl;
        cout << "VictorDoidera_X: " << victorDoidera.getCoords('c')[0] << endl;
        cout << "VictorDoidera_Y: " << victorDoidera.getCoords('c')[1] << endl;
        cout << "VictorDoidera_Z: " << victorDoidera.getCoords('c')[2] << endl;
        cout << pow(2,2) << endl;
        double yDaLoucura = -1 *(sqrt(pow(raioVictoral, 2) - pow(victorDoidera.getCoords('c')[2],2)) - victorDoidera.getCoords('c')[0]);
        cout << "Y da Loucura: " << yDaLoucura << endl;



        Point centroDaLoucura = Point(victor.getCoords('c')[0], yDaLoucura, 0.0, 'c');
        Vector3D sonOfVictor = Vector3D(victor.getCoords('c'), centroDaLoucura.getCoords('c'));
        double yTombado = raioVictoral + yDaLoucura;
        Vector3D sonOfZ = Vector3D( {victor.getCoords('c')[0] , yTombado, 0.0} , centroDaLoucura.getCoords('c'));
        double anguloDaLoucura = sonOfVictor.angle(sonOfZ);
        return anguloDaLoucura;
    };
};

void Molecule::getBonds(){
    string symbol1, symbol2;
    for (int i = 0 ; i < (int) this->molecule.size(); i++){
        for (int j = i; j < (int) this->molecule.size(); j++){
            double length = this->bondLength(i, j);
            symbol1 = this->molecule[i].getAtomicSymbol();
            symbol2 = this->molecule[j].getAtomicSymbol();
            PeriodicTable table = PeriodicTable();
            double radii = this->VDWRatio * (table.getCovalentRadii(symbol1) + table.getCovalentRadii(symbol2));
            if (length <= radii){
                if (i != j){
                    StraightSegment bond = StraightSegment(this->molecule[i].getPoint(), this->molecule[j].getPoint());
                    this->bonds.push_back(pair < vector <int> , StraightSegment > {vector <int> {i,j}, bond});
                };
            };
        };
    };
};

void Molecule::getAngles(){
    for (int i = 0; i < (int) bonds.size(); i++){
        int atom1 = this->bonds[i].first[0];
        int atom2 = this->bonds[i].first[1];
        for (int j = i; j < (int) this->bonds.size(); j++){
            int atom3 = this->bonds[j].first[0];
            int atom4 = this->bonds[j].first[1];
            if (atom1 == atom3){
                if (atom2 != atom4){
                    Angle angle = Angle(this->molecule[atom2].getPoint(), this->molecule[atom1].getPoint(), this->molecule[atom4].getPoint());
                    this->angles.push_back(pair < vector <int>, Angle > {vector <int> {atom2, atom1, atom4}, angle});
                };
            } else if (atom1 == atom4){
                if (atom2 != atom3){
                    Angle angle = Angle(this->molecule[atom2].getPoint(), this->molecule[atom1].getPoint(), this->molecule[atom3].getPoint());
                    this->angles.push_back(pair < vector <int>, Angle > {vector <int> {atom2, atom1, atom3}, angle});
                };
            } else if (atom2 == atom3){
                if (atom1 != atom4){
                    Angle angle = Angle(this->molecule[atom1].getPoint(), this->molecule[atom2].getPoint(), this->molecule[atom4].getPoint());
                    this->angles.push_back(pair < vector <int>, Angle > {vector <int> {atom1, atom2, atom4}, angle});
                };
            } else if (atom2 == atom4){
                if (atom1 != atom3){
                    Angle angle = Angle(this->molecule[atom1].getPoint(), this->molecule[atom2].getPoint(), this->molecule[atom3].getPoint());
                    this->angles.push_back(pair < vector <int>, Angle > {vector <int> {atom1, atom2, atom3}, angle});
                };
            };
        };
    };
};

vector <Atom> Molecule::getMoleculeVector(){
    return this->molecule;
};

void Molecule::getDihedrals(){
    for (int i = 0; i < (int) this->angles.size(); i++){
        int atom1 = this->angles[i].first[0];
        int atom2 = this->angles[i].first[1];
        int atom3 = this->angles[i].first[2];
        for (int j = i; j < (int) this->angles.size(); j++){
            int atom4 = this->angles[j].first[0];
            int atom5 = this->angles[j].first[1];
            int atom6 = this->angles[j].first[2];
            if (atom2 != atom5){
                if (atom2 == atom4 && atom5 == atom3){
                    Torsion dihedral = Torsion(this->molecule[atom1].getPoint(), this->molecule[atom2].getPoint(), this->molecule[atom3].getPoint(), this->molecule[atom6].getPoint());
                    dihedrals.push_back(pair < vector <int>, Torsion > {vector <int> {atom1, atom2, atom3, atom6}, dihedral});
                } else if (atom2 == atom4 && atom5 == atom1){
                    Torsion dihedral = Torsion(this->molecule[atom3].getPoint(), this->molecule[atom2].getPoint(), this->molecule[atom1].getPoint(), this->molecule[atom6].getPoint());
                    dihedrals.push_back(pair < vector <int>, Torsion > {vector <int> {atom3, atom2, atom1, atom6}, dihedral});
                } else if (atom2 == atom6 && atom5 == atom3){
                    Torsion dihedral = Torsion(this->molecule[atom1].getPoint(), this->molecule[atom2].getPoint(), this->molecule[atom3].getPoint(), this->molecule[atom4].getPoint());
                    dihedrals.push_back(pair < vector <int>, Torsion > {vector <int> {atom1, atom2, atom3, atom4}, dihedral});
                } else if (atom2 == atom6 && atom5 == atom1){
                    Torsion dihedral = Torsion(this->molecule[atom3].getPoint(), this->molecule[atom2].getPoint(), this->molecule[atom1].getPoint(), this->molecule[atom4].getPoint());
                    dihedrals.push_back(pair < vector <int>, Torsion > {vector <int> {atom3, atom2, atom1, atom4}, dihedral});
                };
            };
        };
    };
};

vector<double> Molecule::minNmaxValue(vector <double> v){
    vector<double> minMAX(2);
    minMAX.at(0) = v.at(0);
    minMAX.at(1) = v.at(0);
    for(int i = 1; i < (int) v.size(); i++){
        if(v.at(i) < minMAX.at(0)){
            minMAX.at(0) = v.at(i);
        };
        if(v.at(i) > minMAX.at(1)){
            minMAX.at(1) = v.at(i);
        };
    };
    return minMAX;
};

Molecule::Molecule(){
    this->molecule = AtomList();
    this->multiplicity = 1;
    this->VDWRatio = 1.3;
    this->charge = 0;
};

Molecule::~Molecule(){
    molecule.clear();
    chargePoint.clear();
    bonds.clear();
    angles.clear();
    dihedrals.clear();
    multiplicity = 0;
    charge = 0;
};

double Molecule::RMSD(Molecule MOL2)
{
    double rmsd = 0.0;
    for (int i = 0; i < (int) this->molecule.size(); i++)
    {
        rmsd += pow(this->molecule[i].getPos()[0] - MOL2.getAtomObj(i).getPos()[0], 2);
        rmsd += pow(this->molecule[i].getPos()[1] - MOL2.getAtomObj(i).getPos()[1], 2);
        rmsd += pow(this->molecule[i].getPos()[2] - MOL2.getAtomObj(i).getPos()[2], 2);
    }
    return sqrt(rmsd / this->molecule.size());
}

void Molecule::addChargePoints(double xPos, double yPos, double zPos, double charge){
    ChargePoint cp(xPos, yPos, zPos, charge);
    this->chargePoint.push_back(cp);
};

void Molecule::addChargePoints(ChargePoint cp){
    this->chargePoint.push_back(cp);
};

void Molecule::addAtom(string atomSymbol, double xPos, double yPos, double zPos, double atomicCharge, bool freezeCode_){
    Atom atom(atomSymbol, xPos, yPos, zPos, atomicCharge, freezeCode_);
    this->molecule.push_back(atom);
};

void Molecule::addAtom(int atomNumber, double xPos, double yPos, double zPos, double atomicCharge, bool freezeCode_){
    Atom atom(atomNumber, xPos, yPos, zPos, atomicCharge, freezeCode_);
    this->molecule.push_back(atom);
};

void Molecule::addAtom(Atom atom){
    this->molecule.push_back(atom);
};

vector <string> Molecule::getAtom(int number, bool symbol){
    vector<string> atomString(4);
    Atom atom = this->molecule.at(number-1);

    if(symbol == 0){
        atomString.at(0) = to_string(atom.getAtomicNumber());
    }else{
        atomString.at(0) = atom.getAtomicSymbol();
    };

    atomString.at(1) = to_string(atom.getX());
    atomString.at(2) = to_string(atom.getY());
    atomString.at(3) = to_string(atom.getZ());

    return atomString;
};

Atom& Molecule::getAtomObj(int number){
    return this->molecule[number];
}

ChargePoint Molecule::getChargePointsObj(int number){
    return this->chargePoint[number];
}

void Molecule::setVDWRatio(double VDWratio){
    this->VDWRatio = VDWratio;
    this->doIRC();
};

void Molecule::setCharge(int charge){
    this->charge = charge;
};

double Molecule::getCharge(){
    return this->charge;
};

void Molecule::setMultiplicity(int multiplicity){
    this->multiplicity = multiplicity;
};

int Molecule::getMultiplicity(){
    return this->multiplicity;
};

double Molecule::getVDWRatio(){
    return this->VDWRatio;
};

vector< vector<string> > Molecule::getMolecule(bool symbol){
    vector< vector<string> > moleculeString;
    for (int i = 1; i < (int) this->molecule.size()+1; i++){
        vector <string> atom = this->getAtom(i, symbol);
        moleculeString.push_back(atom);
    };
    return moleculeString;
};

vector< vector<string> > Molecule::getChargePoints(){
    vector< vector<string> > cps;
    for (int i=0; i < (int) this->chargePoint.size(); i++){
        vector<string> cp(4);
        cp.at(0) = to_string(this->chargePoint.at(i).getX());
        cp.at(1) = to_string(this->chargePoint.at(i).getY());
        cp.at(2) = to_string(this->chargePoint.at(i).getZ());
        cp.at(3) = to_string(this->chargePoint.at(i).getCharge());
        cps.push_back(cp);
    };
    return cps;
};

long Molecule::getSize(){
    return this->molecule.size();
};

void Molecule::normalizeCPs(int norm){
    for (int i=0; i < (int) this->chargePoint.size(); i++){
        double charge = this->chargePoint.at(i).getCharge();
        this->chargePoint.at(i).setCharge(charge/norm);
    };
};

Point Molecule::getMassCenter(){
    vector <double> massVector;
    vector <double> coordX;
    vector <double> coordY;
    vector <double> coordZ;
    for (int i = 0; i < (int) this->molecule.size(); i++){
        massVector.push_back(this->molecule.at(i).getAtomicMass());
        coordX.push_back(this->molecule.at(i).getX());
        coordY.push_back(this->molecule.at(i).getY());
        coordZ.push_back(this->molecule.at(i).getZ());
    };
    Point temp = MassCenter(massVector, coordX, coordY, coordZ).getMassCenter();
    return temp;
};

void Molecule::spinMolecule(double angle, Vector3D spinVector){
    for (int i = 0; i < (int) this->molecule.size(); i++){
        this->molecule[i].rotationAxis(angle, spinVector);
    };
};

void Molecule::spinMolecule(double angle, char axis){
    if (axis == 'x') {
        Vector3D spinVector = Vector3D({1.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
        this->spinMolecule(angle , spinVector);
    } else if (axis == 'y'){
        Vector3D spinVector = Vector3D({0.0, 1.0, 0.0}, {0.0, 0.0, 0.0});
        this->spinMolecule(angle , spinVector);
    } else {
        Vector3D spinVector = Vector3D({0.0, 0.0, 1.0}, {0.0, 0.0, 0.0});
        this->spinMolecule(angle , spinVector);
    };
};

void Molecule::translation(Vector3D traslationVector){
    for(int i = 0; i < (int) this->molecule.size(); i++){
        this->molecule.at(i).translation(traslationVector);
    };
    if (this->chargePoint.size() != 0){
        for(int i = 0; i < (int) this->chargePoint.size(); i++){
            this->chargePoint.at(i).translation(traslationVector);
        };
    };
};

void Molecule::moveMassCenter(double x, double y, double z){
    Vector3D traslationVector = Vector3D({x, y, z}, this->getMassCenter().getCoords('c'));
    this->translation(traslationVector);
};

void Molecule::moveTail(int atomNumber, double x, double y, double z){
    Vector3D traslationVector = Vector3D({x, y, z}, this->molecule.at(atomNumber).getPos());
    this->translation(traslationVector);
};

Vector3D Molecule::unitVector(Vector3D vector)
{
    double magnitude = vector.magnitude();
    return Vector3D({vector.getVector()[0]/magnitude, vector.getVector()[1]/magnitude, vector.getVector()[2]/magnitude}, {0.0, 0.0, 0.0});
}

double Molecule::bondLength(int atomN1, int atomN2){
    Vector3D bond = Vector3D(this->molecule[atomN1].getPos(), this->molecule[atomN2].getPos());
    return bond.magnitude();
};

double Molecule::valenceAngle(int atomN1, int atomN2, int atomN3){
    Vector3D bond1 = Vector3D(this->molecule[atomN1].getPos(), this->molecule[atomN2].getPos());
    Vector3D bond2 = Vector3D(this->molecule[atomN3].getPos(), this->molecule[atomN2].getPos());
    return bond1.angle(bond2);
};

double Molecule::torsion(int atomN1, int atomN2, int atomN3, int atomN4){
    Vector3D bond1 = Vector3D(this->molecule[atomN2].getPos(), this->molecule[atomN1].getPos());
    Vector3D bond2 = Vector3D(this->molecule[atomN2].getPos(), this->molecule[atomN3].getPos());
    Vector3D bond3 = Vector3D(this->molecule[atomN3].getPos(), this->molecule[atomN4].getPos());
    Vector3D semi_normal1 = bond1.crossProduct(bond2) / sin(bond1.angle(bond2, 'r'));
    Vector3D semi_normal2 = bond3.crossProduct(bond2) / sin(bond3.angle(bond2, 'r'));
    double angleD = semi_normal1.angle(semi_normal2);
    double signal_ = semi_normal1.dotProduct(bond3);
    int signal;
    if (signal_ > 0){
        signal = 1;
    } else {
        signal = -1;
    };
    return signal * angleD;
};

void Molecule::toXYZ(string fileName){

    if (fileName.substr(fileName.find_last_of(".") + 1) != "xyz"){
        fileName = fileName.substr(0, fileName.find_last_of(".")) + ".xyz";
    };

    ofstream file;
    file.open(fileName);
    file << this->molecule.size() << endl;
    file << "XYZ file generated by MoleKing!" << endl;
    for (int i = 0; i < static_cast<int>(this->molecule.size()); i++){
        file << std::left << std::setw(5) << this->molecule[i].getAtomicSymbol() << " "
             << std::fixed << std::setw(12) << this->molecule[i].getX()
             << std::setw(12) << this->molecule[i].getY()
             << std::setw(12) << this->molecule[i].getZ() << std::endl;
    }
    file.close();
}

string toLower(const std::string& str) {
    std::string result = str;
    std::transform(result.begin(), result.end(), result.begin(), ::tolower);
    return result;
}

void Molecule::toGJF(string fileName, string method, string basis, string addKeywords, string midKeywords, string endKeywords, int charge, int multiplicity, bool zmatrix, vector<double> EField)
{    
    auto extension = fileName.substr(fileName.find_last_of(".") + 1);
    if (extension != "gjf" && extension != "com") {
        fileName = fileName.substr(0, fileName.find_last_of(".")) + ".gjf";
    }

    if ( (this->chargePoint.size() == 0) && (midKeywords != ""))
    {
        throw invalid_argument("Please, just use midKeywords if you have charge points, polarizability and electric field.");
    }
    
    if ( (EField.size() > 0) && (toLower(addKeywords).find(toLower("polar")) == string::npos) && (midKeywords != ""))
    {
        throw invalid_argument("Please, just use midKeywords if you have charge points, polarizability and electric field.");
    }

    if (EField.size() > 0) {
        if (!zmatrix) 
        {
            throw invalid_argument("Electric field can only be applied to zmatrix. Set zmatrix=True");
        }

        if (toLower(addKeywords).find(toLower("Field=Read")) == string::npos) 
        {
            addKeywords += " Field=Read";
        }
    }
        
    if (charge != 0){
        this->charge = charge;
    };

    if (multiplicity != 1){
        this->multiplicity = multiplicity;
    };

    if (this->chargePoint.size() != 0){
        if (addKeywords.find("charge") == string::npos) {
            addKeywords = addKeywords + " charge";
        }
    }

    ofstream file;
    file.open(fileName);
    file << "%chk=" << fileName.substr(0, fileName.find_last_of(".")) << ".chk" << endl;
    file << "#p " << method << "/" << basis << " " << addKeywords << endl << endl;
    file << "Gaussian job generated by MoleKing!" << endl << endl;
    file << this->charge << " " << this->multiplicity << endl;

    if (!zmatrix)
    {
        for (int i = 0; i < static_cast<int>(this->molecule.size()); i++){
            file << left << setw(5) << this->molecule[i].getAtomicSymbol() << " "
                << fixed << setw(12) << this->molecule[i].getX()
                << setw(12) << this->molecule[i].getY()
                << setw(12) << this->molecule[i].getZ() << endl;
        }
    }
    
    else 
    {
        this->reorderMolecule();
        this->doIRC();
        int sizeMolecule = this->molecule.size();

        if (sizeMolecule > 0)
        {   
            file << left << setw(5) << this->molecule[0].getAtomicSymbol() << endl;

            if (sizeMolecule > 1)
            {   
                file << left << setw(5) <<  this->molecule[1].getAtomicSymbol() << " "
                << fixed << setw(5) << "1" << " "
                << fixed << setw(15) << to_string(this->bondLength(0, 1)) << endl;

                if (sizeMolecule > 2)
                {      
                    file << left << setw(5) << this->molecule[2].getAtomicSymbol() << " "
                    << fixed << setw(5) << "1" << " "
                    << fixed << setw(15) << to_string(this->bondLength(0, 2)) << " "
                    << fixed << setw(5) << "2" << " "
                    << fixed << setw(15) << to_string(this->valenceAngle(1, 0, 2)) << endl;

                    if (sizeMolecule > 3)
                    {
                        for (int i = 3; (int) i < sizeMolecule; i++)
                        {   
                            file << left << setw(5) << this->molecule[i].getAtomicSymbol() << " "
                            << fixed << setw(5) << to_string(i-2) << " "
                            << fixed << setw(15) << to_string(this->bondLength(i-3,i)) << " "
                            << fixed << setw(5) << to_string(i-1) << " "
                            << fixed << setw(15) << to_string(this->valenceAngle(i, i-3, i-2)) << " "
                            << fixed << setw(5) << to_string(i) << " "
                            << fixed << setw(15) << to_string(this->torsion(i, i-3, i-2,i-1)) << endl;
                        }
                    }
                }
            }
        }
    };
    
    if (this->chargePoint.size() != 0)
    {   
        if ( (EField.size() > 0) && (toLower(addKeywords).find("polar") != string::npos) && (midKeywords != ""))
        {
            file << endl << " " << midKeywords << endl << endl;

            for (int i = 0; i < static_cast<int>(this->chargePoint.size()); i++){
                file << " " <<fixed << setw(12) << this->chargePoint[i].getX()
                    << " " <<fixed << setw(12) << this->chargePoint[i].getY()
                    << " " <<fixed << setw(12) << this->chargePoint[i].getZ() 
                    << " " <<fixed << setw(12) << this->chargePoint[i].getCharge() << endl;
            }

            endKeywords += " " + to_string(EField[0]) + " " + to_string(EField[1]) + " " + to_string(EField[2]);
            file << endl << " " << endKeywords << endl;
        }

        else 
        {
            
            file << endl;

            for (int i = 0; i < static_cast<int>(this->chargePoint.size()); i++)
            {
                file << " " <<fixed << setw(12) << this->chargePoint[i].getX()
                    << " " <<fixed << setw(12) << this->chargePoint[i].getY()
                    << " " <<fixed << setw(12) << this->chargePoint[i].getZ() 
                    << " " <<fixed << setw(12) << this->chargePoint[i].getCharge() << endl;
            }

            if (endKeywords != "")
            {
                if (EField.size() > 0)
                {   
                    endKeywords += " " + to_string(EField[0]) + " " + to_string(EField[1]) + " " + to_string(EField[2]);
                    file << endl << " " << endKeywords << endl;
                }
                else
                {
                    file << endl << " " << endKeywords << endl;
                }
            }
            else
            {
                if (EField.size() > 0)
                {
                    endKeywords += " " + to_string(EField[0]) + " " + to_string(EField[1]) + " " + to_string(EField[2]);
                    file << endl << " " << endKeywords << endl;
                }
                else
                {
                    file << endl;
                }
            }
        }
    }
    else
    {   
        if(EField.size() > 0)
        {
            endKeywords += " " + to_string(EField[0]) + " " + to_string(EField[1]) + " " + to_string(EField[2]);
            file << endl << " " << endKeywords << endl;
        }
        else
        {
            file << endl << " " << endKeywords << endl;
        }
    }

    file << endl;    
    file.close();
}

void Molecule::alignMolecule(char axis){

    Eigen::MatrixXd coords(this->molecule.size(), 3);

    for (size_t i = 0; i < this->molecule.size(); ++i)
    {
        coords.row(i) = Eigen::Vector3d(this->molecule[i].getX(), this->molecule[i].getY(), this->molecule[i].getZ());
    }

    //perform PCA

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(coords.transpose() * coords);
    Eigen::MatrixXd principalAxis = solver.eigenvectors().col(2);  // The principal axis is the eigenvector with the largest eigenvalue
    Eigen::Matrix3d rotationMatrix = Molecule::getRotationMatrix(principalAxis, axis);

    for (size_t i = 0; i < this->molecule.size(); ++i)
    {
        Eigen::Vector3d pos = Eigen::Vector3d(this->molecule[i].getX(), this->molecule[i].getY(), this->molecule[i].getZ());
        pos = rotationMatrix * pos;
        this->molecule[i].setX(pos(0));
        this->molecule[i].setY(pos(1));
        this->molecule[i].setZ(pos(2));
    }

    this->moveMassCenter(0,0,0);
};

Eigen::Matrix3d Molecule::getRotationMatrix(Eigen::Vector3d pAxis, char mkAxis)
{   

    Eigen::Vector3d axis = pAxis.normalized();
    Eigen::Vector3d xAxis(3);

    // if mkAxis == x

    if (mkAxis == 'x')
    {
        xAxis << 1.0, 0.0, 0.0;
    }

    else if (mkAxis == 'y')
    {
        xAxis << 0.0, 1.0, 0.0;
    }

    else if (mkAxis == 'z')
    {
        xAxis << 0.0, 0.0, 1.0;
    }

    else
    {
        throw invalid_argument("Invalid axis. Please, use x, y or z.");
    }

    if (axis.isApprox(xAxis)){
        return Eigen::Matrix3d::Identity();
    }

    Eigen::Vector3d v = axis.cross(xAxis);

    double s = v.norm();
    double c = axis.dot(xAxis);

    Eigen::Matrix3d vx;
    vx << 0, -v(2), v(1),
          v(2), 0, -v(0),
          -v(1), v(0), 0;

    Eigen::Matrix3d rotationMatrix = Eigen::Matrix3d::Identity() + vx + vx * vx * (1 - c) / (s * s);
    return rotationMatrix;
};



void Molecule::reorderMolecule()
{
    //create a function that reorders the molecule based on the distance with the first atom

    vector <int> temp;
    vector <double> distances;
    vector <double> atom1 = this->molecule[0].getPos();
    for (int i = 1; i < (int) this->molecule.size(); i++){
        vector <double> atom2 = this->molecule[i].getPos();
        distances.push_back(Vector3D(atom1, atom2).magnitude());
    }
    vector <double> distancesCopy = distances;
    sort(distances.begin(), distances.end());
    for (int i = 0; i < (int) distances.size(); i++){
        for (int j = 0; j < (int) distancesCopy.size(); j++){
            if (distances[i] == distancesCopy[j]){
                temp.push_back(j+1);
                distancesCopy[j] = -1;
                break;
            }
        }
    }
    AtomList tempMolecule;
    tempMolecule.push_back(this->molecule[0]);
    for (int i = 0; i < (int) temp.size(); i++){
        tempMolecule.push_back(this->molecule[temp[i]]);
    }
    this->molecule = tempMolecule;

};

// Define the method doZMatrix for the Molecule class
void Molecule::doIRC(){
    this->bonds.clear();
    this->angles.clear();
    this->dihedrals.clear();
    this->getBonds();
    this->getAngles();
    this->getDihedrals();
};


vector < vector <int> > Molecule::getIRCBonds(){
    if (this->bonds.size() == 0){
        this->doIRC();
    }    
    vector < vector <int> > temp;
    for (int i = 0; i < (int) this->bonds.size(); i++){
        temp.push_back(this->bonds[i].first);
    }
    return temp;
};

vector < vector <int> > Molecule::getIRCAngles(){
    if (this->angles.size() == 0){
        this->doIRC();
    }
    vector < vector <int> > temp;
    for (int i = 0; i < (int) this->angles.size(); i++){
        temp.push_back(this->angles[i].first);
    }
    return temp;
};

vector < vector <int> > Molecule::getIRCDihedrals(){
    if (this->dihedrals.size() == 0){
        this->doIRC();
    }
    //vector < vector <int> > temp(this->dihedrals.size());
    vector < vector <int> > temp;
    for (int i = 0; i < (int) this->dihedrals.size(); i++){
        temp.push_back(this->dihedrals[i].first);
    }
    return temp;
};

vector <Atom> Molecule::moleculeList(){
    return this->molecule;
};

Molecule::iterator Molecule::begin(){
    return this->molecule.begin();
};

Molecule::iterator Molecule::end(){
    return this->molecule.end();
};

void Molecule::removeAtom(int atomNumber){
    this->molecule.erase(this->molecule.begin() + atomNumber);
};

void Molecule::clear()
{
    this->molecule.clear();
};

void Molecule::removeAtom(Atom atom){
    for (int i = 0; i < (int) this->molecule.size(); i++){
        if (atom == this->molecule[i]){
            this->molecule.erase(this->molecule.begin() + i);
            break;
        };
    };
};

string Molecule::toStr(){
    vector <pair <string, int> > s;
    string result = "Molecule ";
    s.push_back(pair <string, int> {this->molecule[0].getAtomicSymbol(), 1});
    for (int i = 1; i < (int) this->molecule.size(); i++){
        string symbol = this->molecule[i].getAtomicSymbol();
        for (int j = 0; j < (int) s.size(); j++){
            if (symbol == s[j].first){
                int value = s[j].second + 1;
                s.at(j) = (pair<string, int> {symbol, value});
                break;
            } else if (j == (int) s.size()-1) {
                s.push_back(pair <string, int> {symbol, 1});
                break;
            } else {
                continue;
            };
        };
    };
    for (int i = 0; i < (int) s.size(); i++){
        result = result + s[i].first + "_{"  + to_string(s[i].second) + "}";
    };
    result = result + ", with charge " + to_string(this->charge);
    if (this->multiplicity != 0){
        result = result + " and multiplicity " + to_string(this->multiplicity) ;
    };
    if (this->chargePoint.size() != 0){
        result = result + " and with " + to_string(this->chargePoint.size()) + " charge points";
    };
    return result;
};

bool Molecule::operator==(Molecule mol){
    if ((int) this->molecule.size() ==(int) mol.getSize()){
        for (int i = 0; i < (int) this->molecule.size(); i++){
            if (this->molecule[i] == mol.getAtomObj(i)){
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

bool Molecule::operator!=(Molecule mol){
    if ((int) this->molecule.size() == (int) mol.getSize()){
        for (int i = 0; i < (int) this->molecule.size(); i++){
            if (this->molecule[i] == mol.getAtomObj(i)){
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

bool Molecule::operator<(Molecule mol){
    if ((int) this->molecule.size() < (int) mol.getSize()){
        return 1;
    } else {
        return 0;
    };
};

bool Molecule::operator>(Molecule mol){
    if ((int) this->molecule.size() > (int) mol.getSize()){
        return 1;
    } else {
        return 0;
    };
};

void Molecule::removeElement(string element){

    Molecule temp = Molecule();

    for (int i = 0; i < this->molecule.size(); i++){
        if (this->molecule[i].getAtomicSymbol() != element){
            temp.addAtom(this->molecule[i]);
        };
    }

    this->molecule = temp.molecule;
};

Molecule Molecule::copy(){
    Molecule temp = Molecule();
    for (int i = 0; i < (int) this->molecule.size(); i++){
        temp.addAtom(this->molecule[i]);
    };
    for (int i = 0; i < (int) this->chargePoint.size(); i++){
        temp.addChargePoints(this->chargePoint[i]);
    };
    return temp;
};

double Molecule::getMolecularMass(){
    double molMass = 0;
    for (int i = 0; i < (int) this->molecule.size(); i++) {
        molMass = molMass + this->molecule[i].getAtomicMass();
    }
    return molMass;
};
