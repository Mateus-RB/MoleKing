//   MoleKing //
//
//   File:        [SampleMolecules.hpp]
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


#ifndef SampleMolecules_hpp
#define SampleMolecules_hpp

#include <stdio.h>
#include <map>
#include <iterator>
#include <string.h>
#include <vector>
#include <sstream>
#include <iostream>

using namespace std;

class SampleMolecules{
    private:
    map<string, vector<vector<string>>> sampleMoleculesMap;
    //map<string, int> sampleMoleculesSize;
    vector<string> customSplit(string str, char separator = ' ');
    vector<vector<string>> makeMol(vector<char> myChar, char sep);

    public:
    SampleMolecules();
    ~SampleMolecules();
    vector<vector<string>> getSampleMoleculesVector(string mol);
};

#endif /* SampleMolecules_hpp */
