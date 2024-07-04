//   MoleKing //
//
//   File:        [main.cpp]
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

#include <iostream>
#include <string>
#include <math.h>
#include "../include/Eigen/Eigenvalues"
#include "MoleKing.hpp"
#include "myMath/MassCenter.hpp"
#include "chemicalUnits/AtomicScale.hpp"
#include "myMath/Geometry.hpp"
#include "myMath/Matrix.hpp"
#include "chemicalUnits/Molecule.hpp"
#include "berny/Hessian.hpp"
#include "chemicalUnits/SupraMolecule.hpp"
#include "chemicalUnits/OPLSff.hpp"
#include "myMath/Vectors.hpp"
#include "outputProcess/G16Process.hpp"
#include "outputProcess/Psi4Process.hpp"
#include "chemicalUnits/pov.cpp"
#include "chemicalUnits/pov.hpp"
using namespace std;

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/pytypes.h>
namespace py = pybind11;
PYBIND11_MODULE(MoleKing, m) {
    
    py::class_<PeriodicTable>(m, "PeriodicTable", "This class creates a virtual Periodic Table.")
        .def(py::init())
        .def("getAtomicNumber", &PeriodicTable::getAtomicNumber)
        .def("getAtomicMass", &PeriodicTable::getAtomicMass)
        .def("getSymbol", &PeriodicTable::getSymbol)
        .def("getCovalentRadii", &PeriodicTable::getCovalentRadii)
        .def("getColor", &PeriodicTable::getColor);
    
    py::class_<Atom>(m, "Atom", "This class creates a atom variable type allowing for the usage in python like a primitive type.")
        .def(py::init<int, double, double, double, double, bool>(), py::arg("atomicNumber"), py::arg("xPos"), py::arg("yPos"), py::arg("zPos"), py::arg("atomicCharge") = 0.0, py::arg("freezeCode_") = 0)
        .def(py::init<string, double, double, double, double, bool>(), py::arg("atomicSymbol"), py::arg("xPos"), py::arg("yPos"), py::arg("zPos"), py::arg("atomicCharge") = 0.0, py::arg("freezeCode_") = 0)
        .def("getAtomicMass", &Atom::getAtomicMass)
        .def("getAtomicSymbol", &Atom::getAtomicSymbol)
        .def("getAtomicNumber", &Atom::getAtomicNumber)
        .def("getAtomicRadio", &Atom::getAtomicRadio)
        .def("getX", &Atom::getX)
        .def("getY", &Atom::getY)
        .def("getZ", &Atom::getZ)
        .def("getAtomicCharge", &Atom::getAtomicCharge)
        .def("setX", &Atom::setX)
        .def("setY", &Atom::setY)
        .def("setZ", &Atom::setZ)
        .def("setAtomicCharge", &Atom::setAtomicCharge)
        .def("__eq__", &Atom::operator==)
        .def("__ne__", &Atom::operator==)
        .def("__lt__", &Atom::operator<)
        .def("__gt__", &Atom::operator>)
        .def("__le__", &Atom::operator<=)
        .def("__ge__", &Atom::operator>=)
        .def("__str__", &Atom::toStr)
        .def("setNewPos", &Atom::setNewPos)
        .def("getPos", &Atom::getPos)
        .def("translation", &Atom::translation)
        .def("rotationAxis", &Atom::rotationAxis);
    
    
    py::class_<ChargePoint>(m, "ChargePoint", "This class creates a charge point variable type allowing for the usage in python like a primitive type.")
        .def(py::init<double, double, double, double>(), py::arg("xPos"), py::arg("yPos"), py::arg("zPos"), py::arg("charge"))
        .def("getX", &ChargePoint::getX)
        .def("getY", &ChargePoint::getY)
        .def("getZ", &ChargePoint::getZ)
        .def("setX", &ChargePoint::setX)
        .def("setY", &ChargePoint::setY)
        .def("setZ", &ChargePoint::setZ)
        .def("getCharge", &ChargePoint::getCharge)
        .def("setCharge", &ChargePoint::setCharge)
        .def("setNewPos", &ChargePoint::setNewPos)
        .def("getPos", &ChargePoint::getPos)
        .def("__eq__", &ChargePoint::operator==)
        .def("__ne__", &ChargePoint::operator!=)
        .def("__lt__", &ChargePoint::operator<)
        .def("__gt__", &ChargePoint::operator>)
        .def("__le__", &ChargePoint::operator<=)
        .def("__ge__", &ChargePoint::operator>=)
        .def("__str__", &ChargePoint::toStr)
        .def("translation", &ChargePoint::translation)
        .def("rotationAxis", &ChargePoint::rotationAxis);
    
    py::class_<Molecule>(m, "Molecule", "This class creates a molecule variable type allowing for the usage in python like a primitive type.")
        .def(py::init())
        .def("addChargePoints", (void (Molecule::*)(double, double, double, double)) &Molecule::addChargePoints, "This method add a charge point in a existent molecule.")
        .def("addChargePoints", (void (Molecule::*)(ChargePoint)) &Molecule::addChargePoints, "This method add a charge point in a existent molecule.")
        .def("addAtom", (void (Molecule::*)(string, double, double, double, double, bool)) &Molecule::addAtom, py::arg("atomSymbol"), py::arg("xPos"), py::arg("yPos"), py::arg("zPos"),py::arg("atomicCharge")=0.0, py::arg("freezeCode_")=0)
        .def("addAtom", (void (Molecule::*)(int, double, double, double, double, bool)) &Molecule::addAtom, py::arg("atomNumber"), py::arg("xPos"), py::arg("yPos"), py::arg("zPos"), py::arg("atomicCharge")=0.0, py::arg("freezeCode_")=0)
        .def("addAtom", (void (Molecule::*)(Atom)) &Molecule::addAtom)
        .def("removeAtom", (void (Molecule::*)(int)) &Molecule::removeAtom)
        .def("removeAtom", (void (Molecule::*)(Atom)) &Molecule::removeAtom)
        .def("getAtom", &Molecule::getAtom, py::arg("number")=0, py::arg("symbol")=0)
        .def("__getitem__", &Molecule::getAtomObj, py::return_value_policy::reference)
        .def("__str__", &Molecule::toStr)
        .def("setCharge", &Molecule::setCharge)
        .def("getCharge", &Molecule::getCharge)
        .def("__eq__", &Molecule::operator==)
        .def("__ne__", &Molecule::operator!=)
        .def("__lt__", &Molecule::operator<)
        .def("__gt__", &Molecule::operator>)
        .def("__len__", &Molecule::getSize)
        .def("__iter__", [](Molecule &mol) {return py::make_iterator(mol.begin(), mol.end());}, py::keep_alive<0, 1>())
        .def("setMultiplicity", &Molecule::setMultiplicity)
        .def("setVDWRatio", &Molecule::setVDWRatio)
        .def("getVDWRatio", &Molecule::getVDWRatio)
        .def("getMultiplicity", &Molecule::getMultiplicity)
        .def("copy", &Molecule::copy)
        .def("clear", &Molecule::clear)
        .def("normChargePoints", &Molecule::normalizeCPs)
        .def("getMolecule", &Molecule::getMolecule)
        .def("getChargePoints", &Molecule::getChargePoints)
        .def("getMassCenter", &Molecule::getMassCenter)
        .def("spinMolecule", (void (Molecule::*)(double, Vector3D)) &Molecule::spinMolecule)
        .def("spinMolecule", (void (Molecule::*)(double, char)) &Molecule::spinMolecule)
        .def("translation", &Molecule::translation)
        .def("moveMassCenter", &Molecule::moveMassCenter, py::arg("x")=0, py::arg("y")=0, py::arg("z")=0)
        .def("moveTail", &Molecule::moveTail, py::arg("atomNumber"), py::arg("x")=0, py::arg("y")=0, py::arg("z")=0)
        .def("bondLength", &Molecule::bondLength)
        .def("valenceAngle", &Molecule::valenceAngle)
        .def("torsion", &Molecule::torsion)
        .def("getIRCBonds", &Molecule::getIRCBonds)
        .def("getIRCAngles", &Molecule::getIRCAngles)
        .def("getIRCDihedrals", &Molecule::getIRCDihedrals)
        .def("removeElement", &Molecule::removeElement)
        .def("toXYZ", &Molecule::toXYZ, py::arg("fileName") = "MK_Molecule.xyz")
        .def("RMSD", &Molecule::RMSD, py::arg("MOL2"))
        .def("toGJF", &Molecule::toGJF, py::arg("fileName") = "MK_Molecule.gjf", py::arg("method") = "B3LYP", py::arg("basis") = "6-311g(d)", py::arg("addKeywords") = "", py::arg("midKeywords") = "", py::arg("endKeywords") = "", py::arg("charge") = 0, py::arg("multiplicity") = 1, py::arg("zmatrix") = 0, py::arg("EField") = vector<double> {})
        .def("alignMolecule", &Molecule::alignMolecule, py::arg("axis") = 'x')
        .def("getMM", &Molecule::getMolecularMass);

    py::class_<SupraMolecule>(m, "SupraMolecule", "This class creates a set of molecules variable type allowing for the usage in python like a primitive type.")
        .def(py::init())
        .def("addMolecule", &SupraMolecule::addMolecule)
        .def("addAtomToMolecule", (void (SupraMolecule::*)(int, Atom)) &SupraMolecule::addAtomToMolecule)
        .def("addAtomToMolecule", (void (SupraMolecule::*)(int, string, double, double, double)) &SupraMolecule::addAtomToMolecule)
        .def("__str__", &SupraMolecule::toStr)
        .def("getMolecule", &SupraMolecule::getMolecule)
        .def("setMultiplicity", &SupraMolecule::setMultiplicity)
        .def("getMultiplicity", &SupraMolecule::getMultiplicity)
        .def("getCharge", &SupraMolecule::getCharge)
        .def("getSize", &SupraMolecule::getSize)
        .def("getMassCenter", &SupraMolecule::getMassCenter)
        .def("translation", &SupraMolecule::translation)
        .def("moveMassCenter", &SupraMolecule::moveMassCenter)
        .def("moveTail", &SupraMolecule::moveTail)
        .def("spinSupraMolecule", (void (SupraMolecule::*)(double, char)) &SupraMolecule::spinSupraMolecule)
        .def("spinSupraMolecule", (void (SupraMolecule::*)(double, Vector3D)) &SupraMolecule::spinSupraMolecule)
        .def("getIRCBonds", &SupraMolecule::getIRCBonds)
        .def("getIRCAngles", &SupraMolecule::getIRCAngles)
        .def("getIRCDihedrals", &SupraMolecule::getIRCDihedrals)
        .def("removeAtom", (void (SupraMolecule::*)(int, int)) &SupraMolecule::removeAtom)
        .def("removeAtom", (void (SupraMolecule::*)(Atom)) &SupraMolecule::removeAtom)
        .def("removeElement", (void (SupraMolecule::*)(int, string)) &SupraMolecule::removeElement)
        .def("removeElement", (void (SupraMolecule::*)(string)) &SupraMolecule::removeElement)
        .def("removeMolecule", (void (SupraMolecule::*)(int)) &SupraMolecule::removeMolecule)
        .def("removeMolecule", (void (SupraMolecule::*)(Molecule)) &SupraMolecule::removeMolecule)
        .def("__iter__", [](SupraMolecule &sMol) {return py::make_iterator(sMol.begin(), sMol.end());}, py::keep_alive<0, 1>())
        .def("__eq__", &SupraMolecule::operator==)
        .def("__ne__", &SupraMolecule::operator!=)
        .def("__lt__", &SupraMolecule::operator<)
        .def("__gt__", &SupraMolecule::operator>)
        .def("getSupraMM", &SupraMolecule::getSupraMolecularMass);;
    
    py::class_<Point>(m, "Point", "This class creates a point variable type allowing for the usage in python like a primitive type.")
        .def(py::init())
        .def(py::init<double, double, double, char>(), py::arg("coord1"), py::arg("coord2"), py::arg("coord3"), py::arg("typeCoord") = 'c')
        .def("__eq__", &Point::operator==)
        .def("getCoords", &Point::getCoords)
        .def("setCoord", &Point::setCoord)
        .def("setCoords", &Point::setCoords, py::arg("newValues"), py::arg("typeCoord") = 'c')
        .def("translation", &Point::translation)
        .def("rotation3D", &Point::rotationVector)
        .def("__str__", &Point::toStr, py::arg("spaceType") = 'c');
    
    py::class_<SphericalCoords>(m, "SphericalCoords", "This class allows the interchange between Cartesian and spherical coordinates.")
        .def(py::init<double, double, double, char>(), py::arg("coord1"), py::arg("coord2"), py::arg("coord3"), py::arg("spaceType") = 'c')
        .def("toCartesian", &SphericalCoords::toCartesian)
        .def("toSpherical", &SphericalCoords::toSpherical);
    
    py::class_<Vector3D>(m, "Vector3D", "This class creates a vector (xi + yj + zk) variable type allowing for the usage in python like a primitive type.")
        .def(py::init< vector<double>, vector<double>>(), py::arg("pointA"), py::arg("pointB") = vector <double> {0.0, 0.0, 0.0})
        .def("__abs__", &Vector3D::magnitude)
        .def("getVector", &Vector3D::getVector)
        .def("normalize", &Vector3D::normalize)
        .def("__invert__", &Vector3D::conjugate)
        .def("__div__", &Vector3D::operator/)
        .def("__mul__", &Vector3D::operator*)
        .def("__add__", &Vector3D::operator+)
        .def("__sub__", &Vector3D::operator-)
        .def("__str__", &Vector3D::toStr)
        .def("crossProduct", &Vector3D::crossProduct)
        .def("dotProduct", &Vector3D::dotProduct)
        .def("angle", &Vector3D::angle, py::arg("vectorB"), py::arg("unit") = 'd')
        .def("unitVectorValue", &Vector3D::axisValue)
        .def("magnitude", &Vector3D::magnitude);
    
    py::class_<Matrix>(m, "Matrix", "This class creates a Matrix variable type allowing for the usage in python like a primitive type.")
        .def(py::init< vector < vector<double> > >())
        .def(py::init< int, int >())
        .def(py::init())
        .def("setMatrix", &Matrix::setMatrix)
        .def("add", &Matrix::sum)
        .def("multip", (Matrix (Matrix::*)(double)) &Matrix::multiplication)
        .def("multip", (Matrix (Matrix::*)(Matrix)) &Matrix::multiplication)
        .def("determinant", &Matrix::determinant)
        .def("replace", &Matrix::replace)
        .def("__getitem__", &Matrix::getLine)
        .def("elem", &Matrix::element)
        .def("show", &Matrix::print)
        .def("__str__", &Matrix::toStr);

    py::class_<PovRay>(m, "PovRay", "This class creates a PovRay file.")
        .def(py::init<const Molecule&>(), py::arg("mol"));
        
    py::class_<G16LOGfile>(m, "G16LOGfile", "This class is experimental and under development.")
        .def(py::init< string, bool, bool, int>(), py::arg("filePath"), py::arg("polarAsw") = 0, py::arg("tdAsw") = 0, py::arg("link") = 0)
        .def("getDate", &G16LOGfile::getDate)
        .def("getEnergy", &G16LOGfile::getEnergy)
        .def("getBasis", &G16LOGfile::getBasis)
        .def("getMethod", &G16LOGfile::getMethod)
        .def("getMolecule", &G16LOGfile::getMolecule)
        .def("getOrbitals", &G16LOGfile::getOrbitals)
        .def("getTransitions", &G16LOGfile::getTransitions, py::arg("index") = 0)
        .def("getDipole", &G16LOGfile::getDipole, py::arg("axis") = "tot")
        .def("getFrequency", &G16LOGfile::getFrequency)
        .def("getNLO", &G16LOGfile::getNLO, py::arg("Orientation") = "input")
        .def("getAlpha", &G16LOGfile::getAlpha,py::arg("orientation") = "Dipole", py::arg("unit") = "esu", py::arg("frequency") = 0)
        .def("getBeta", &G16LOGfile::getBeta,py::arg("orientation") = "Dipole", py::arg("unit") = "esu", py::arg("frequency") = 0, py::arg("BSHG") = 0)
        .def("getGamma", &G16LOGfile::getGamma,py::arg("orientation") = "Dipole", py::arg("unit") = "esu", py::arg("frequency") = 0, py::arg("GSHG") = 0)
        .def("getHOMO", [](G16LOGfile &self, int index) {
            auto values = self.getHOMO(index);
            if (values.size() == 1)
            {
                return py::cast(values[0]);
            }
            return py::cast(values);
        }, py::arg("index") = 0)
        .def("getLUMO", [](G16LOGfile &self, int index) {
            auto values = self.getLUMO(index);
            if (values.size() == 1)
            {
                return py::cast(values[0]);
            }
            return py::cast(values);
        }, py::arg("index") = 0)
        .def("__str__", &G16LOGfile::toStr);   

    py::class_<Psi4OUTfile>(m, "Psi4OUTfile", "This class is experimental and under development.")
        .def(py::init< string >(), py::arg("filePath"))
        .def("getMolecule", &Psi4OUTfile::getMolecule)
        .def("getMul", &Psi4OUTfile::getMul)
        .def("getCharge", &Psi4OUTfile::getCharge)
        .def("__str__", &Psi4OUTfile::toStr);    
};
