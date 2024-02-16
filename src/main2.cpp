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
#include <../include/Eigen/Eigenvalues>
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

using namespace std;

int main(int argc, char* argv[]){
    /*
    */
      if (argc < 2)
      {
            cout << "MoleKing" << " Version " << MoleKing_VERSION_MAJOR << "." << MoleKing_VERSION_MINOR << "." << MoleKing_VERSION_PATCH << endl;
            cout << "Usage: " << argv[0] << " parameter" << endl;
            return 1;
      } 
  
      else if ((string) argv[1] == "water")
      {
            Atom A1("O", 0.3326888, 0.5289030, 0.0000);
            Atom A2(1, 1.2926887, 0.5291441, 0.0000);
            Atom A3(1, 0.0124615, 1.43391, 0.0000);
            Molecule water = Molecule();
            water.addAtom(A1);
            water.addAtom(A2);
            water.addAtom(A3);
            cout << water.toStr() << endl;
            vector < vector <int> > b = water.getIRCBonds();
            for (int i = 0; i < b.size(); i++)
            {
                  for (int j = 0; j < b[i].size(); j++)
                  {
                        cout << b[i][j] << endl;
                  }
            }
        
            return 0;
      } 
  
      else if ((string) argv[1] == "--version")
      {
            cout << "MoleKing" << " Version " << MoleKing_VERSION_MAJOR << "." << MoleKing_VERSION_MINOR << "." << MoleKing_VERSION_PATCH << endl;
            return 0;
      } 
  
      else if ((string) argv[1] == "methanol")
      {
            Atom A1("C",-3.89100000, 0.48400000, 0.00000000);
            Atom A2("H",-3.50500000, 0.96700000, 0.87200000);
            Atom A3("H",-4.96100000, 0.51800000, 0.01900000);
            Atom A4("H",-3.56700000,-0.53600000,-0.01600000);
            Atom A5("O",-3.41400000, 1.15800000,-1.16800000);
            Atom A6("H",-3.76000000, 0.72300000,-1.95100000);
            Molecule met = Molecule();
            met.addAtom(A1);
            met.addAtom(A2);
            met.addAtom(A3);
            met.addAtom(A4);
            met.addAtom(A5);
            met.addAtom(A6);
            vector < vector <int> > b = met.getIRCAngles();
            for ( auto i = 0; i < b.size(); i++ )
            {
                  cout << b[i][0] << " " << b[i][1] << " " << b[i][2] << endl;
            };
            return 0;
      }   

      else if ((string) argv[1] == "charges")
      {
            Atom A1("C",-3.89100000, 0.48400000, 0.00000000);            
            cout << "Charge of atom B4 -> " << A1.getAtomicCharge() << endl;
            A1.setAtomicCharge(1.0);
            cout << "Charge of atom After -> " << A1.getAtomicCharge() << endl;

            cout << "\n" << endl;

            Molecule M1 = Molecule();
            M1.addAtom("H", -3.50500000, 0.96700000, 0.87200000);
            cout << "Charge of atom B4 -> " << M1.getAtomObj(0).getAtomicCharge() << endl;
            M1.getAtomObj(0).setAtomicCharge(1.0);
            cout << "Charge of atom After -> " << M1.getAtomObj(0).getAtomicCharge() << endl;
      }
 
      else{
            int AN1, AN2, AN3;
            double X1, X2, X3, Y1, Y2, Y3, Z1, Z2, Z3;
            cout <<"Add 3 Atoms" << endl;
            cout << "----------------" << endl;
            cout << "Atomic number of first Atom" << endl;
            cin >> AN1;
            cout << "X coordinated of first Atom" << endl;
            cin >> X1;
            cout << "Y coordinated of first Atom" << endl;
            cin >> Y1;
            cout << "Z coordinated of first Atom" << endl;
            cin >> Z1;
            cout << "----------------" << endl;
            cout << "Atomic number of second Atom" << endl;
            cin >> AN2;
            cout << "X coordinated of second Atom" << endl;
            cin >> X2;
            cout << "Y coordinated of second Atom" << endl;
            cin >> Y2;
            cout << "Z coordinated of second Atom" << endl;
            cin >> Z2;
            cout << "----------------" << endl;
            cout << "Atomic number of third Atom" << endl;
            cin >> AN3;
            cout << "X coordinated of third Atom" << endl;
            cin >> X3;
            cout << "Y coordinated of third Atom" << endl;
            cin >> Y3;
            cout << "Z coordinated of third Atom" << endl;
            cin >> Z3;
            Atom A1(AN1, X1, Y1, Z1);
            Atom A2(AN2, X2, Y2, Z2);
            Atom A3(AN3, X3, Y3, Z3);
            Molecule M = Molecule();
            M.addAtom(A1);
            M.addAtom(A2);
            M.addAtom(A3);
            cout << M.getVDWRatio() << endl;
            M.setVDWRatio(2.0);
            cout << M.getVDWRatio() << endl;
            cout << "The molecule is:" << endl;
            cout << M.toStr() << endl;
            return 0;
      };
};


