//   MoleKing //
//
//   File:        [Molecule.cpp]
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
    } else {
        Point victor = Point(cart[0], cart[1], cart[2], 'c');
        Point victorDoidera = Point(cart[0], cart[1], cart[2], 'c');
        Vector3D yAxis = Vector3D({0.0, 1.0, 0.0}, {0.0, 0.0, 0.0});
        victorDoidera.rotationVector(180, yAxis);
        Vector3D victores = Vector3D(victor.getCoords('c'), victorDoidera.getCoords('c'));
        double linhaDaLoucura = victores.magnitude();
        double raioVictoral = linhaDaLoucura/2;
        double yDaLoucura = -1 *(sqrt(pow(raioVictoral, 2) - pow(victorDoidera.getCoords('c')[2],2)) - victorDoidera.getCoords('c')[1]);
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

Atom Molecule::getAtomObj(int number){
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

void Molecule::standardOrientation(){
    vector<int> biggerDistance = this->molecularAxis();
    this->moveTail(biggerDistance[0]);
    double angle1 = this->angleToSpinInAref(biggerDistance[1], 'y');
    this->spinMolecule(angle1, 'y');
    vector <double> zspin = this->molecule[biggerDistance[1]].getPos();
    vector <double> zspinSpherical = SphericalCoords(zspin[0], zspin[1], zspin[2], 'c').toSpherical();
    this->spinMolecule(zspinSpherical[2], 'z');
    double angle2 = this->angleToSpinInAref(biggerDistance[1], 'x');
    this->spinMolecule(angle2, 'x');
    this->spinMolecule(90, 'z');
    this->spinMolecule(90, 'y');
    this->moveMassCenter();
};

vector <double> Molecule::standardOrientationPath(){
    vector<int> biggerDistance = this->molecularAxis();
    this->moveTail(biggerDistance[0]);
    double angle1 = this->angleToSpinInAref(biggerDistance[1], 'y');
    this->spinMolecule(angle1, 'y');
    vector <double> zspin = this->molecule[biggerDistance[1]].getPos();
    vector <double> zspinSpherical = SphericalCoords(zspin[0], zspin[1], zspin[2], 'c').toSpherical();
    this->spinMolecule(zspinSpherical[2], 'z');
    double angle2 = this->angleToSpinInAref(biggerDistance[1], 'x');
    this->spinMolecule(angle2, 'x');
    this->spinMolecule(90, 'z');
    this->spinMolecule(90, 'y');
    this->moveMassCenter();
    return vector <double> {angle1, zspinSpherical[2], angle2};
};

Eigen::Matrix<double, 3, 3> Molecule::moleculeTensor(){
    Point centerOfMass = this->getMassCenter();
    double Ixx = 0.0;
    double Iyy = 0.0;
    double Izz = 0.0;
    double Ixy = 0.0;
    double Ixz = 0.0;
    double Iyz = 0.0;

    //double atomMass = this->molecule.at(0).getAtomicMass();

    for (int i = 0; i < (int) this->molecule.size(); i++){
        double atomMass = this->molecule.at(i).getAtomicMass();        
        double dx = this->molecule.at(i).getX() - centerOfMass.getCoords('c')[0];
        double dy = this->molecule.at(i).getY() - centerOfMass.getCoords('c')[1];
        double dz = this->molecule.at(i).getZ() - centerOfMass.getCoords('c')[2];

        Ixx += atomMass * (dy*dy + dz*dz);
        Iyy += atomMass * (dx*dx + dz*dz);
        Izz += atomMass * (dx*dx + dy*dy);
        Ixy += atomMass * dx*dy * -1;
        Ixz += atomMass * dx*dz * -1;
        Iyz += atomMass * dy*dz * -1;
    };

    Matrix inertiaTensor = Matrix({{Ixx,Ixy,Ixz},
                                   {Ixy,Iyy,Iyz},
                                   {Ixz,Iyz,Izz} });

    
    Eigen::Matrix<double, 3, 3> A;

    Eigen::Matrix<double, 1, 3> Position;

    //Position << -0.280418, 0.000000,-0.000026;

    A << Ixx, Ixy, Ixz,
         Ixy, Iyy, Iyz,
         Ixz, Iyz, Izz;

    Eigen::EigenSolver<Eigen::Matrix<double, 3, 3> > s(A); 
    Eigen::Matrix<double, 3, 3> evecs = s.eigenvectors().real(); //! Eigenvectors are the columns of evecs.
    
    //cout << "Position: " << endl << Position << endl;
    //cout << "\n" << endl;
    //cout << "Position Transpose: " << endl << Position.transpose() << endl;
    //cout << "\n" << endl;
    //cout << "Eigenvalues:" << endl << s.eigenvalues() << endl;
    //cout << "\n" << endl;
    //cout << "Eigenvectors: " << endl << evecs << endl;
    //cout << "\n" << endl;
    //cout << "Eigenvectors Tranpose" << endl << evecs.transpose() << endl;
    //cout << "\n" << endl;
    //cout << "M_Multi w/out Tranpose" << endl << Position*evecs << endl;
    //cout << "\n" << endl;

    //cout << "\n\n";

    return evecs;
}

void Molecule::stdOrientation(){
    Eigen::Matrix<double, 3, 3> evecs = this->moleculeTensor();

    if (float(evecs.determinant()) == float(-1)){
        evecs(0,2) = evecs(0,2) * -1;
        evecs(1,2) = evecs(1,2) * -1;
        evecs(2,2) = evecs(2,2) * -1;
    }

    else{
        cout << "Error: could not make a rotation matrix while adopting the standard orientation" << endl;
    };   
    
    this->moveMassCenter(0,0,0);

    for (int i = 0; i < (int) this->molecule.size(); i++ ){
        Eigen::Matrix<double, 3, 1> temp;
        temp << 0.0, 
                0.0, 
                0.0;

        temp(0,0) = this->molecule.at(i).getX();
        temp(1,0) = this->molecule.at(i).getY();
        temp(2,0) = this->molecule.at(i).getZ();

        //cout << "\nTemp \n\n";
        //cout << temp;      
        //cout << "\nTemp Transpose\n\n";
        //cout << temp.transpose();
        //cout << "\nEvecs Transpose\n\n";
        //cout<< evecs.transpose();
        //cout << "\n";

        Eigen::Matrix<double, 1, 3> newPosition = (temp.transpose() * evecs.transpose()).transpose();

        float atomX = newPosition(0,0);
        float atomY = newPosition(0,1);
        float atomZ = newPosition(0,2);

        this->molecule.at(i).setX(atomX);
        this->molecule.at(i).setY(atomY);
        this->molecule.at(i).setZ(atomZ);        
//
        temp(0,0) = 0;
        temp(1,0) = 0;
        temp(2,0) = 0;
    }

    return;
}


vector <int> Molecule::molecularAxis(){
    int j = 0;
    vector <int> temp;
    temp.push_back(0);
    temp.push_back(0);
    double distance = 0;
    while(j < (int) this->molecule.size()){
        for(int i = j+1; i < (int) this->molecule.size(); i++){
            vector<double> atomCoord1 = this->molecule.at(j).getPos();
            vector<double> atomCoord2 = this->molecule.at(i).getPos();
            double dist = Vector3D(atomCoord1, atomCoord2).magnitude();
            if(dist > distance){
                distance = dist;
                temp.at(0) = j;
                temp.at(1) = i;
            };
        };
        j++;
    };
    return temp;
};

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

void Molecule::toGJF(string fileName, string method, string basis, string addKeywords,  string endKeywords, int charge, int multiplicity){
    
    if (fileName.substr(fileName.find_last_of(".") + 1) != "gjf" && fileName.substr(fileName.find_last_of(".") + 1) != "com"){
        fileName = fileName.substr(0, fileName.find_last_of(".")) + ".gjf";
    };
    
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
    for (int i = 0; i < static_cast<int>(this->molecule.size()); i++){
        file << left << setw(5) << this->molecule[i].getAtomicSymbol() << " "
             << fixed << setw(12) << this->molecule[i].getX()
             << setw(12) << this->molecule[i].getY()
             << setw(12) << this->molecule[i].getZ() << endl;
    }
    
    if (this->chargePoint.size() != 0){
        if (endKeywords != ""){
            file << endl << " " << endKeywords << endl << endl;
        }
        else {
            file << endl;
        }
        for (int i = 0; i < static_cast<int>(this->chargePoint.size()); i++){
            file << " " <<fixed << setw(12) << this->chargePoint[i].getX()
                 << " " <<fixed << setw(12) << this->chargePoint[i].getY()
                 << " " <<fixed << setw(12) << this->chargePoint[i].getZ() 
                 << " " <<fixed << setw(12) << this->chargePoint[i].getCharge() << endl;
        }
    }

    else {
        file << endl << " " << endKeywords << endl;
    }

    file << endl;    
    file.close();
}

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

void Molecule::removeElement(string element){
    for (int i = 0; i < (int) this->molecule.size(); i++) {
        if (this->molecule[i].getAtomicSymbol() == element){
            this->molecule.erase(this->molecule.begin() + i);
        };
    };
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
