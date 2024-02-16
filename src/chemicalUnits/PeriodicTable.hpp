//   MoleKing //
//
//   File:        [PeriodicTable.hpp]
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


#ifndef PeriodicTable_hpp
#define PeriodicTable_hpp

#include <stdio.h>
#include <map>
#include <iterator>
#include <string>

using namespace std;


class PeriodicTable{
    private:
    map<string, int> symbolMap;
    map<string, double> massMap;
    map<string, double> radiiMap;
    map<string, string> ColorMap;
    public:
    PeriodicTable();
    ~PeriodicTable();
    int getAtomicNumber(string symbol);
    double getAtomicMass(string symbol);
    string getSymbol(int atomicNumber);
    double getCovalentRadii(string symbol);
    string getColor(string symbol);
};

#endif /* PeriodicTable_hpp */
