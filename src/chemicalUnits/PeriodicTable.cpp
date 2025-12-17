//   MoleKing //
//
//   File:        [PeriodicTable.cpp]
//
//   Author(s):   ['LEEDMOL Reserch Group']
//   Site(s):     ['https://www.researchgate.net/lab/LEEDMOL-Heibbe-Cristhians-Research-Group-Heibbe-Cristhian-B-De-Oliveira']
//   Email(s):    ['heibbe@ufg.br']
//   Credits:     ['Copyright © 2023 LEEDMOL. All rights reserved.']
//   Date:        ['17.01.2023']
//   Version:     ['1.6.0']
//   Status:      ['Development']
//   Language:    ['C++','Python']
//   Description: ['A python module written in C++ for theoretical chemestry']


#include "PeriodicTable.hpp"

PeriodicTable::PeriodicTable(){
    //Atomic number table
    this->symbolMap.insert(pair<string, int>("H", 1));
    this->symbolMap.insert(pair<string, int>("He", 2));
    this->symbolMap.insert(pair<string, int>("Li", 3));
    this->symbolMap.insert(pair<string, int>("Be", 4));
    this->symbolMap.insert(pair<string, int>("B", 5));
    this->symbolMap.insert(pair<string, int>("C", 6));
    this->symbolMap.insert(pair<string, int>("N", 7));
    this->symbolMap.insert(pair<string, int>("O", 8));
    this->symbolMap.insert(pair<string, int>("F", 9));
    this->symbolMap.insert(pair<string, int>("Ne", 10));
    this->symbolMap.insert(pair<string, int>("Na", 11));
    this->symbolMap.insert(pair<string, int>("Mg", 12));
    this->symbolMap.insert(pair<string, int>("Al", 13));
    this->symbolMap.insert(pair<string, int>("Si", 14));
    this->symbolMap.insert(pair<string, int>("P", 15));
    this->symbolMap.insert(pair<string, int>("S", 16));
    this->symbolMap.insert(pair<string, int>("Cl", 17));
    this->symbolMap.insert(pair<string, int>("Ar", 18));
    this->symbolMap.insert(pair<string, int>("K", 19));
    this->symbolMap.insert(pair<string, int>("Ca", 20));
    this->symbolMap.insert(pair<string, int>("Sc", 21));
    this->symbolMap.insert(pair<string, int>("Ti", 22));
    this->symbolMap.insert(pair<string, int>("V", 23));
    this->symbolMap.insert(pair<string, int>("Cr", 24));
    this->symbolMap.insert(pair<string, int>("Mn", 25));
    this->symbolMap.insert(pair<string, int>("Fe", 26));
    this->symbolMap.insert(pair<string, int>("Co", 27));
    this->symbolMap.insert(pair<string, int>("Ni", 28));
    this->symbolMap.insert(pair<string, int>("Cu", 29));
    this->symbolMap.insert(pair<string, int>("Zn", 30));
    this->symbolMap.insert(pair<string, int>("Ga", 31));
    this->symbolMap.insert(pair<string, int>("Ge", 32));
    this->symbolMap.insert(pair<string, int>("As", 33));
    this->symbolMap.insert(pair<string, int>("Se", 34));
    this->symbolMap.insert(pair<string, int>("Br", 35));
    this->symbolMap.insert(pair<string, int>("Kr", 36));
    this->symbolMap.insert(pair<string, int>("Rb", 37));
    this->symbolMap.insert(pair<string, int>("Sr", 38));
    this->symbolMap.insert(pair<string, int>("Y", 39));
    this->symbolMap.insert(pair<string, int>("Zr", 40));
    this->symbolMap.insert(pair<string, int>("Nb", 41));
    this->symbolMap.insert(pair<string, int>("Mo", 42));
    this->symbolMap.insert(pair<string, int>("Tc", 43));
    this->symbolMap.insert(pair<string, int>("Ru", 44));
    this->symbolMap.insert(pair<string, int>("Rh", 45));
    this->symbolMap.insert(pair<string, int>("Pd", 46));
    this->symbolMap.insert(pair<string, int>("Ag", 47));
    this->symbolMap.insert(pair<string, int>("Cd", 48));
    this->symbolMap.insert(pair<string, int>("In", 49));
    this->symbolMap.insert(pair<string, int>("Sn", 50));
    this->symbolMap.insert(pair<string, int>("Sb", 51));
    this->symbolMap.insert(pair<string, int>("Te", 52));
    this->symbolMap.insert(pair<string, int>("I", 53));
    this->symbolMap.insert(pair<string, int>("Xe", 54));
    this->symbolMap.insert(pair<string, int>("Cs", 55));
    this->symbolMap.insert(pair<string, int>("Ba", 56));
    this->symbolMap.insert(pair<string, int>("La", 57));
    this->symbolMap.insert(pair<string, int>("Ce", 58));
    this->symbolMap.insert(pair<string, int>("Pr", 59));
    this->symbolMap.insert(pair<string, int>("Nd", 60));
    this->symbolMap.insert(pair<string, int>("Pm", 61));
    this->symbolMap.insert(pair<string, int>("Sm", 62));
    this->symbolMap.insert(pair<string, int>("Eu", 63));
    this->symbolMap.insert(pair<string, int>("Gd", 64));
    this->symbolMap.insert(pair<string, int>("Tb", 65));
    this->symbolMap.insert(pair<string, int>("Dy", 66));
    this->symbolMap.insert(pair<string, int>("Ho", 67));
    this->symbolMap.insert(pair<string, int>("Er", 68));
    this->symbolMap.insert(pair<string, int>("Tm", 69));
    this->symbolMap.insert(pair<string, int>("Yb", 70));
    this->symbolMap.insert(pair<string, int>("Lu", 71));
    this->symbolMap.insert(pair<string, int>("Hf", 72));
    this->symbolMap.insert(pair<string, int>("Ta", 73));
    this->symbolMap.insert(pair<string, int>("W", 74));
    this->symbolMap.insert(pair<string, int>("Re", 75));
    this->symbolMap.insert(pair<string, int>("Os", 76));
    this->symbolMap.insert(pair<string, int>("Ir", 77));
    this->symbolMap.insert(pair<string, int>("Pt", 78));
    this->symbolMap.insert(pair<string, int>("Au", 79));
    this->symbolMap.insert(pair<string, int>("Hg", 80));
    this->symbolMap.insert(pair<string, int>("Tl", 81));
    this->symbolMap.insert(pair<string, int>("Pb", 82));
    this->symbolMap.insert(pair<string, int>("Bi", 83));
    this->symbolMap.insert(pair<string, int>("Po", 84));
    this->symbolMap.insert(pair<string, int>("At", 85));
    this->symbolMap.insert(pair<string, int>("Rn", 86));
    this->symbolMap.insert(pair<string, int>("Fr", 87));
    this->symbolMap.insert(pair<string, int>("Ra",88));
    this->symbolMap.insert(pair<string, int>("Ac", 89));
    this->symbolMap.insert(pair<string, int>("Th", 90));
    this->symbolMap.insert(pair<string, int>("Pa", 91));
    this->symbolMap.insert(pair<string, int>("U",  92));
    this->symbolMap.insert(pair<string, int>("Np", 93));
    this->symbolMap.insert(pair<string, int>("Pu", 94));
    this->symbolMap.insert(pair<string, int>("Am", 95));
    this->symbolMap.insert(pair<string, int>("Cm", 96));
    this->symbolMap.insert(pair<string, int>("Bk", 97));
    this->symbolMap.insert(pair<string, int>("Cf", 98));
    this->symbolMap.insert(pair<string, int>("Es", 99));
    this->symbolMap.insert(pair<string, int>("Fm", 100));
    this->symbolMap.insert(pair<string, int>("Md", 101));
    this->symbolMap.insert(pair<string, int>("No", 102));
    this->symbolMap.insert(pair<string, int>("Lr", 103));
    this->symbolMap.insert(pair<string, int>("Rf", 104));
    this->symbolMap.insert(pair<string, int>("Db", 105));
    this->symbolMap.insert(pair<string, int>("Sg", 106));
    this->symbolMap.insert(pair<string, int>("Bh", 107));
    this->symbolMap.insert(pair<string, int>("Hs", 108));
    this->symbolMap.insert(pair<string, int>("Mt", 109));
    this->symbolMap.insert(pair<string, int>("Ds", 110));
    this->symbolMap.insert(pair<string, int>("Rg", 111));
    this->symbolMap.insert(pair<string, int>("Cn", 112));
    this->symbolMap.insert(pair<string, int>("Nh", 113));
    this->symbolMap.insert(pair<string, int>("Fl", 114));
    this->symbolMap.insert(pair<string, int>("Mc", 115));
    this->symbolMap.insert(pair<string, int>("Lv", 116));
    this->symbolMap.insert(pair<string, int>("Og", 117));
    this->symbolMap.insert(pair<string, int>("Ts", 118));
    this->symbolMap.insert(pair<string, int>("00", 0));
    this->symbolMap.insert(pair<string, int>("dice_ghost_label", 0));

    // Mass Table, All masses were taken from: https://doi.org/10.1002/rcm.8864. Atoms that were not listed there were taken elsewhere (references in comments aside).
    this->massMap.insert(pair<string, double>("H",  1.007975000));
    this->massMap.insert(pair<string, double>("He", 4.002602000));
    this->massMap.insert(pair<string, double>("Li", 6.967500000));
    this->massMap.insert(pair<string, double>("Be", 9.012183100));
    this->massMap.insert(pair<string, double>("B",  10.813500000));
    this->massMap.insert(pair<string, double>("C",  12.010600000));
    this->massMap.insert(pair<string, double>("N",  14.006855000));
    this->massMap.insert(pair<string, double>("O",  15.999400000));
    this->massMap.insert(pair<string, double>("F",  18.998403162));
    this->massMap.insert(pair<string, double>("Ne", 20.179700000));
    this->massMap.insert(pair<string, double>("Na", 22.989769280));
    this->massMap.insert(pair<string, double>("Mg", 24.305500000));
    this->massMap.insert(pair<string, double>("Al", 26.981538400));
    this->massMap.insert(pair<string, double>("Si", 28.085000000));
    this->massMap.insert(pair<string, double>("P",  30.973761998));
    this->massMap.insert(pair<string, double>("S",  32.067500000));
    this->massMap.insert(pair<string, double>("Cl", 35.451500000));
    this->massMap.insert(pair<string, double>("Ar", 39.862500000));
    this->massMap.insert(pair<string, double>("K",  39.098300000));
    this->massMap.insert(pair<string, double>("Ca", 40.078000000));
    this->massMap.insert(pair<string, double>("Sc", 44.955907000));
    this->massMap.insert(pair<string, double>("Ti", 47.867000000));
    this->massMap.insert(pair<string, double>("V",  50.941500000));
    this->massMap.insert(pair<string, double>("Cr", 51.996100000));
    this->massMap.insert(pair<string, double>("Mn", 54.938043000));
    this->massMap.insert(pair<string, double>("Fe", 55.845000000));
    this->massMap.insert(pair<string, double>("Co", 58.933194000));
    this->massMap.insert(pair<string, double>("Ni", 58.693400000));
    this->massMap.insert(pair<string, double>("Cu", 63.546000000));
    this->massMap.insert(pair<string, double>("Zn", 65.380000000));
    this->massMap.insert(pair<string, double>("Ga", 69.723000000));
    this->massMap.insert(pair<string, double>("Ge", 72.630000000));
    this->massMap.insert(pair<string, double>("As", 74.921595000));
    this->massMap.insert(pair<string, double>("Se", 79.971000000));
    this->massMap.insert(pair<string, double>("Br", 79.904000000));
    this->massMap.insert(pair<string, double>("Kr", 83.798000000));
    this->massMap.insert(pair<string, double>("Rb", 85.467800000));
    this->massMap.insert(pair<string, double>("Sr", 87.620000000));
    this->massMap.insert(pair<string, double>("Y",  88.905838000));
    this->massMap.insert(pair<string, double>("Zr", 91.224000000));
    this->massMap.insert(pair<string, double>("Nb", 92.906370000));
    this->massMap.insert(pair<string, double>("Mo", 95.950000000));
    this->massMap.insert(pair<string, double>("Tc", 97.907216000)); // https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
    this->massMap.insert(pair<string, double>("Ru", 101.070000000));
    this->massMap.insert(pair<string, double>("Rh", 102.905490000));
    this->massMap.insert(pair<string, double>("Pd", 106.420000000));
    this->massMap.insert(pair<string, double>("Ag", 107.868200000));
    this->massMap.insert(pair<string, double>("Cd", 112.414000000));
    this->massMap.insert(pair<string, double>("In", 114.818000000));
    this->massMap.insert(pair<string, double>("Sn", 118.710000000));
    this->massMap.insert(pair<string, double>("Sb", 121.760000000));
    this->massMap.insert(pair<string, double>("Te", 127.600000000));
    this->massMap.insert(pair<string, double>("I",  126.904470000));
    this->massMap.insert(pair<string, double>("Xe", 131.293000000));
    this->massMap.insert(pair<string, double>("Cs", 132.905451960));
    this->massMap.insert(pair<string, double>("Ba", 137.327000000));
    this->massMap.insert(pair<string, double>("La", 138.905470000));
    this->massMap.insert(pair<string, double>("Ce", 140.116000000));
    this->massMap.insert(pair<string, double>("Pr", 140.907660000));
    this->massMap.insert(pair<string, double>("Nd", 144.242000000));
    this->massMap.insert(pair<string, double>("Pm", 144.912744000)); // https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
    this->massMap.insert(pair<string, double>("Sm", 150.360000000));
    this->massMap.insert(pair<string, double>("Eu", 151.964000000));
    this->massMap.insert(pair<string, double>("Gd", 157.250000000));
    this->massMap.insert(pair<string, double>("Tb", 158.925354000));
    this->massMap.insert(pair<string, double>("Dy", 162.500000000));
    this->massMap.insert(pair<string, double>("Ho", 164.930329000));
    this->massMap.insert(pair<string, double>("Er", 167.259000000));
    this->massMap.insert(pair<string, double>("Tm", 168.934219000));
    this->massMap.insert(pair<string, double>("Yb", 173.045000000));
    this->massMap.insert(pair<string, double>("Lu", 174.966800000));
    this->massMap.insert(pair<string, double>("Hf", 178.486000000));
    this->massMap.insert(pair<string, double>("Ta", 180.947880000));
    this->massMap.insert(pair<string, double>("W",  183.840000000));
    this->massMap.insert(pair<string, double>("Re", 186.207000000));
    this->massMap.insert(pair<string, double>("Os", 190.230000000));
    this->massMap.insert(pair<string, double>("Ir", 192.217000000));
    this->massMap.insert(pair<string, double>("Pt", 195.084000000));
    this->massMap.insert(pair<string, double>("Au", 196.966570000));
    this->massMap.insert(pair<string, double>("Hg", 200.592000000));
    this->massMap.insert(pair<string, double>("Tl", 204.383500000));
    this->massMap.insert(pair<string, double>("Pb", 207.040000000));
    this->massMap.insert(pair<string, double>("Bi", 208.980400000));
    this->massMap.insert(pair<string, double>("Po", 208.982416000)); // https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
    this->massMap.insert(pair<string, double>("At", 209.987131000)); // https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
    this->massMap.insert(pair<string, double>("Rn", 222.017570000));  // https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
    this->massMap.insert(pair<string, double>("Fr", 223.019731000)); // https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
    this->massMap.insert(pair<string, double>("Ra", 226.025403000)); // https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
    this->massMap.insert(pair<string, double>("Ac", 227.027747000)); // https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
    this->massMap.insert(pair<string, double>("Th", 232.037700000)); 
    this->massMap.insert(pair<string, double>("Pa", 231.035880000));
    this->massMap.insert(pair<string, double>("U",  238.028910000));
    this->massMap.insert(pair<string, double>("Np", 237.048167000)); // https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf 
    this->massMap.insert(pair<string, double>("Pu", 244.064198000)); // https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
    this->massMap.insert(pair<string, double>("Am", 243.061373000)); // https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
    this->massMap.insert(pair<string, double>("Cm", 247.070347000)); // https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
    this->massMap.insert(pair<string, double>("Bk", 247.070299000)); // https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
    this->massMap.insert(pair<string, double>("Cf", 251.079580000)); // https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
    this->massMap.insert(pair<string, double>("Es", 252.082972000)); // https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
    this->massMap.insert(pair<string, double>("Fm", 257.095099000)); // https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
    this->massMap.insert(pair<string, double>("Md", 258.098425000)); // https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
    this->massMap.insert(pair<string, double>("No", 259.101024000)); // https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
    this->massMap.insert(pair<string, double>("Lr", 262.109692000)); // https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
    this->massMap.insert(pair<string, double>("Rf", 263.118313000)); // https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
    this->massMap.insert(pair<string, double>("Db", 262.011437000)); // https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
    this->massMap.insert(pair<string, double>("Sg", 266.012238000)); // https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
    this->massMap.insert(pair<string, double>("Bh", 264.012496000)); // https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
    this->massMap.insert(pair<string, double>("Hs", 269.001341000)); // https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
    this->massMap.insert(pair<string, double>("Mt", 268.001388000)); // https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
    this->massMap.insert(pair<string, double>("Ds", 272.001463000)); // https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
    this->massMap.insert(pair<string, double>("Rg", 272.001535000)); // https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
    this->massMap.insert(pair<string, double>("Cn", 285.177120000)); // https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl
    this->massMap.insert(pair<string, double>("Nh", 284.178730000)); // https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl
    this->massMap.insert(pair<string, double>("Fl", 289.190420000)); // https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl
    this->massMap.insert(pair<string, double>("Mc", 288.192740000)); // https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl
    this->massMap.insert(pair<string, double>("Lv", 293.204490000)); // https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl
    this->massMap.insert(pair<string, double>("Ts", 292.207460000)); // https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl
    this->massMap.insert(pair<string, double>("Og", 294.213920000)); // https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl
    this->massMap.insert(pair<string, double>("00", 0.000000000));
    this->massMap.insert(pair<string, double>("dice_ghost_label", 0.000000000));

    //Covalent Raddi as in https://pubs.acs.org/doi/pdf/10.1021/jp5065819 R1
    this->radiiMap.insert(pair<string, double>("H", 0.320));
    this->radiiMap.insert(pair<string, double>("He", 0.460));
    this->radiiMap.insert(pair<string, double>("Li", 1.330));
    this->radiiMap.insert(pair<string, double>("Be", 1.020));
    this->radiiMap.insert(pair<string, double>("B", 0.850));
    this->radiiMap.insert(pair<string, double>("C", 0.750));
    this->radiiMap.insert(pair<string, double>("N", 0.710));
    this->radiiMap.insert(pair<string, double>("O", 0.630));
    this->radiiMap.insert(pair<string, double>("F", 0.640));
    this->radiiMap.insert(pair<string, double>("Ne", 0.670));
    this->radiiMap.insert(pair<string, double>("Na", 1.550));
    this->radiiMap.insert(pair<string, double>("Mg", 1.390));
    this->radiiMap.insert(pair<string, double>("Al", 1.260));
    this->radiiMap.insert(pair<string, double>("Si", 1.160));
    this->radiiMap.insert(pair<string, double>("P", 1.110));
    this->radiiMap.insert(pair<string, double>("S", 1.030));
    this->radiiMap.insert(pair<string, double>("Cl", 0.990));
    this->radiiMap.insert(pair<string, double>("Ar", 0.960));
    this->radiiMap.insert(pair<string, double>("K", 1.960));
    this->radiiMap.insert(pair<string, double>("Ca", 1.710));
    this->radiiMap.insert(pair<string, double>("Sc", 1.480));
    this->radiiMap.insert(pair<string, double>("Ti", 1.360));
    this->radiiMap.insert(pair<string, double>("V", 1.340));
    this->radiiMap.insert(pair<string, double>("Cr", 1.220));
    this->radiiMap.insert(pair<string, double>("Mn", 1.190));
    this->radiiMap.insert(pair<string, double>("Fe", 1.160));
    this->radiiMap.insert(pair<string, double>("Co", 1.110));
    this->radiiMap.insert(pair<string, double>("Ni", 1.100));
    this->radiiMap.insert(pair<string, double>("Cu", 1.120));
    this->radiiMap.insert(pair<string, double>("Zn", 1.180));
    this->radiiMap.insert(pair<string, double>("Ga", 1.240));
    this->radiiMap.insert(pair<string, double>("Ge", 1.210));
    this->radiiMap.insert(pair<string, double>("As", 1.210));
    this->radiiMap.insert(pair<string, double>("Se", 1.160));
    this->radiiMap.insert(pair<string, double>("Br", 1.140));
    this->radiiMap.insert(pair<string, double>("Kr", 1.170));
    this->radiiMap.insert(pair<string, double>("Rb", 2.100));
    this->radiiMap.insert(pair<string, double>("Sr", 1.850));
    this->radiiMap.insert(pair<string, double>("Y", 1.630));
    this->radiiMap.insert(pair<string, double>("Zr", 1.540));
    this->radiiMap.insert(pair<string, double>("Nb", 1.470));
    this->radiiMap.insert(pair<string, double>("Mo", 1.380));
    this->radiiMap.insert(pair<string, double>("Tc", 1.280));
    this->radiiMap.insert(pair<string, double>("Ru", 1.250));
    this->radiiMap.insert(pair<string, double>("Rh", 1.250));
    this->radiiMap.insert(pair<string, double>("Pd", 1.200));
    this->radiiMap.insert(pair<string, double>("Ag", 1.280));
    this->radiiMap.insert(pair<string, double>("Cd", 1.360));
    this->radiiMap.insert(pair<string, double>("In", 1.420));
    this->radiiMap.insert(pair<string, double>("Sn", 1.400));
    this->radiiMap.insert(pair<string, double>("Sb", 1.400));
    this->radiiMap.insert(pair<string, double>("Te", 1.360));
    this->radiiMap.insert(pair<string, double>("I", 1.330));
    this->radiiMap.insert(pair<string, double>("Xe", 1.310));
    this->radiiMap.insert(pair<string, double>("Cs", 2.320));
    this->radiiMap.insert(pair<string, double>("Ba", 1.960));
    this->radiiMap.insert(pair<string, double>("La", 1.800));
    this->radiiMap.insert(pair<string, double>("Ce", 1.630));
    this->radiiMap.insert(pair<string, double>("Pr", 1.760));
    this->radiiMap.insert(pair<string, double>("Nd", 1.740));
    this->radiiMap.insert(pair<string, double>("Pm", 1.730));
    this->radiiMap.insert(pair<string, double>("Sm", 1.720));
    this->radiiMap.insert(pair<string, double>("Eu", 1.680));
    this->radiiMap.insert(pair<string, double>("Gd", 1.690));
    this->radiiMap.insert(pair<string, double>("Tb", 1.680));
    this->radiiMap.insert(pair<string, double>("Dy", 1.670));
    this->radiiMap.insert(pair<string, double>("Ho", 1.660));
    this->radiiMap.insert(pair<string, double>("Er", 1.650));
    this->radiiMap.insert(pair<string, double>("Tm", 1.640));
    this->radiiMap.insert(pair<string, double>("Yb", 1.700));
    this->radiiMap.insert(pair<string, double>("Lu", 1.620));
    this->radiiMap.insert(pair<string, double>("Hf", 1.520));
    this->radiiMap.insert(pair<string, double>("Ta", 1.460));
    this->radiiMap.insert(pair<string, double>("W", 1.370));
    this->radiiMap.insert(pair<string, double>("Re", 1.310));
    this->radiiMap.insert(pair<string, double>("Os", 1.290));
    this->radiiMap.insert(pair<string, double>("Ir", 1.220));
    this->radiiMap.insert(pair<string, double>("Pt", 1.230));
    this->radiiMap.insert(pair<string, double>("Au", 1.240));
    this->radiiMap.insert(pair<string, double>("Hg", 1.330));
    this->radiiMap.insert(pair<string, double>("Tl", 1.440));
    this->radiiMap.insert(pair<string, double>("Pb", 1.440));
    this->radiiMap.insert(pair<string, double>("Bi", 1.510));
    this->radiiMap.insert(pair<string, double>("Po", 1.450));
    this->radiiMap.insert(pair<string, double>("At", 1.470));
    this->radiiMap.insert(pair<string, double>("Rn", 1.420));
    this->radiiMap.insert(pair<string, double>("Fr", 2.230));
    this->radiiMap.insert(pair<string, double>("Ra", 2.010));
    this->radiiMap.insert(pair<string, double>("Ac", 1.860));
    this->radiiMap.insert(pair<string, double>("Th", 1.750));
    this->radiiMap.insert(pair<string, double>("Pa", 1.690));
    this->radiiMap.insert(pair<string, double>("U",  1.700));
    this->radiiMap.insert(pair<string, double>("Np", 1.710));
    this->radiiMap.insert(pair<string, double>("Pu", 1.720));
    this->radiiMap.insert(pair<string, double>("Am", 1.660));
    this->radiiMap.insert(pair<string, double>("Cm", 1.660));
    this->radiiMap.insert(pair<string, double>("Bk", 1.680));
    this->radiiMap.insert(pair<string, double>("Cf", 1.680));
    this->radiiMap.insert(pair<string, double>("Es", 1.650));
    this->radiiMap.insert(pair<string, double>("Fm", 1.670));
    this->radiiMap.insert(pair<string, double>("Md", 1.730));
    this->radiiMap.insert(pair<string, double>("No", 1.760));
    this->radiiMap.insert(pair<string, double>("Lr", 1.610));
    this->radiiMap.insert(pair<string, double>("Rf", 1.570));
    this->radiiMap.insert(pair<string, double>("Db", 1.490));
    this->radiiMap.insert(pair<string, double>("Sg", 1.430));
    this->radiiMap.insert(pair<string, double>("Bh", 1.410));
    this->radiiMap.insert(pair<string, double>("Hs", 1.340));
    this->radiiMap.insert(pair<string, double>("Mt", 1.290));
    this->radiiMap.insert(pair<string, double>("Ds", 1.280));
    this->radiiMap.insert(pair<string, double>("Rg", 1.210));
    this->radiiMap.insert(pair<string, double>("Cn", 1.220));
    this->radiiMap.insert(pair<string, double>("Nh", 1.360));
    this->radiiMap.insert(pair<string, double>("Fl", 1.430));
    this->radiiMap.insert(pair<string, double>("Mc", 1.620));
    this->radiiMap.insert(pair<string, double>("Lv", 1.750));
    this->radiiMap.insert(pair<string, double>("Og", 1.650));
    this->radiiMap.insert(pair<string, double>("Ts", 1.570));
    this->radiiMap.insert(pair<string, double>("00", 0));
    this->radiiMap.insert(pair<string, double>("dice_ghost_label", 0));

    // Color Map in HEX format
    this->ColorMap.insert(pair<string, string>("H",         "FFFFFF"));
    this->ColorMap.insert(pair<string, string>("He",        "D9FFFF"));
    this->ColorMap.insert(pair<string, string>("Li",        "CC80FF"));
    this->ColorMap.insert(pair<string, string>("Be",        "C2FF00"));
    this->ColorMap.insert(pair<string, string>("B",         "FFB5B5"));
    this->ColorMap.insert(pair<string, string>("C",         "7C7C7C"));
    this->ColorMap.insert(pair<string, string>("N",         "3050F8"));
    this->ColorMap.insert(pair<string, string>("O",         "FF0D0D"));
    this->ColorMap.insert(pair<string, string>("F",         "90E050"));
    this->ColorMap.insert(pair<string, string>("Ne",        "B3E3F5"));
    this->ColorMap.insert(pair<string, string>("Na",        "AB5CF2"));
    this->ColorMap.insert(pair<string, string>("Mg",        "8AFF00"));
    this->ColorMap.insert(pair<string, string>("Al",        "BFA6A6"));
    this->ColorMap.insert(pair<string, string>("Si",        "F0C8A0"));
    this->ColorMap.insert(pair<string, string>("P",         "FF8000"));
    this->ColorMap.insert(pair<string, string>("S",         "FFFF30"));
    this->ColorMap.insert(pair<string, string>("Cl",        "1FF01F"));
    this->ColorMap.insert(pair<string, string>("Br",        "32ff00"));
    this->ColorMap.insert(pair<string, string>("Ar",        "80D1E3"));
    this->ColorMap.insert(pair<string, string>("K",         "8F40D4"));
    this->ColorMap.insert(pair<string, string>("Ca",        "3DFF00"));
    this->ColorMap.insert(pair<string, string>("Sc",        "E6E6E6"));
    this->ColorMap.insert(pair<string, string>("Ti",        "BFC2C7"));
    this->ColorMap.insert(pair<string, string>("V",         "A6A6AB"));
    this->ColorMap.insert(pair<string, string>("Cr",        "8A99C7"));
    this->ColorMap.insert(pair<string, string>("Mn",        "9C7AC7"));
    this->ColorMap.insert(pair<string, string>("Fe",        "E06633"));
    this->ColorMap.insert(pair<string, string>("Co",        "F090A0"));
    this->ColorMap.insert(pair<string, string>("Ni",        "50D050"));
    this->ColorMap.insert(pair<string, string>("Cu",        "C88033"));
    this->ColorMap.insert(pair<string, string>("Zn",        "7D80B0"));
    this->ColorMap.insert(pair<string, string>("Y",         "94FFFF"));
    this->ColorMap.insert(pair<string, string>("Zr",        "94E0E0"));
    this->ColorMap.insert(pair<string, string>("Nb",        "73C2C9"));
    this->ColorMap.insert(pair<string, string>("Mo",        "54B5B5"));
    this->ColorMap.insert(pair<string, string>("Tc",        "3B9E9E"));
    this->ColorMap.insert(pair<string, string>("Ru",        "248F8F"));
    this->ColorMap.insert(pair<string, string>("Rh",        "0A7D8C"));
    this->ColorMap.insert(pair<string, string>("Pd",        "006985"));
    this->ColorMap.insert(pair<string, string>("Ag",        "C0C0C0"));
    this->ColorMap.insert(pair<string, string>("Cd",        "FFD98F"));
    this->ColorMap.insert(pair<string, string>("In",        "A67573"));
    this->ColorMap.insert(pair<string, string>("Sn",        "668080"));
    this->ColorMap.insert(pair<string, string>("Sb",        "9E63B5"));
    this->ColorMap.insert(pair<string, string>("Te",        "D47A00"));
    this->ColorMap.insert(pair<string, string>("I",         "940094"));
    this->ColorMap.insert(pair<string, string>("Xe",        "429EB0"));
    this->ColorMap.insert(pair<string, string>("Cs",        "57178F"));
    this->ColorMap.insert(pair<string, string>("Ba",        "00C900"));
    this->ColorMap.insert(pair<string, string>("La",        "70D4FF"));
    this->ColorMap.insert(pair<string, string>("Ce",        "FFFFC7"));
    this->ColorMap.insert(pair<string, string>("Pr",        "D9FFC7"));
    this->ColorMap.insert(pair<string, string>("Nd",        "C7FFC7"));
    this->ColorMap.insert(pair<string, string>("Pm",        "A3FFC7"));
    this->ColorMap.insert(pair<string, string>("Sm",        "8FFFC7"));
    this->ColorMap.insert(pair<string, string>("Eu",        "61FFC7"));
    this->ColorMap.insert(pair<string, string>("Gd",        "45FFC7"));
    this->ColorMap.insert(pair<string, string>("Tb",        "30FFC7"));
    this->ColorMap.insert(pair<string, string>("Dy",        "1FFFC7"));
    this->ColorMap.insert(pair<string, string>("Ho",        "00FF9C"));
    this->ColorMap.insert(pair<string, string>("Er",        "00E675"));
    this->ColorMap.insert(pair<string, string>("Tm",        "00D452"));
    this->ColorMap.insert(pair<string, string>("Yb",        "00BF38"));
    this->ColorMap.insert(pair<string, string>("Lu",        "00AB24"));
    this->ColorMap.insert(pair<string, string>("Hf",        "4DC2FF"));
    this->ColorMap.insert(pair<string, string>("Ta",        "4DA6FF"));
    this->ColorMap.insert(pair<string, string>("W",         "2194D6"));
    this->ColorMap.insert(pair<string, string>("Re",        "267DAB"));
    this->ColorMap.insert(pair<string, string>("Os",        "266696"));
    this->ColorMap.insert(pair<string, string>("Ir",        "175487"));
    this->ColorMap.insert(pair<string, string>("Pt",        "D0D0E0"));
    this->ColorMap.insert(pair<string, string>("Au",        "FFD123"));
    this->ColorMap.insert(pair<string, string>("Hg",        "B8B8D0"));
    this->ColorMap.insert(pair<string, string>("Tl",        "A6544D"));
    this->ColorMap.insert(pair<string, string>("Pb",        "575961"));
    this->ColorMap.insert(pair<string, string>("Bi",        "9E4FB5"));
    this->ColorMap.insert(pair<string, string>("Po",        "AB5C00"));
    this->ColorMap.insert(pair<string, string>("At",        "754F45"));
    this->ColorMap.insert(pair<string, string>("Rn",        "428296"));
    this->ColorMap.insert(pair<string, string>("Fr",        "420066"));
    this->ColorMap.insert(pair<string, string>("Ra",        "007D00"));
    this->ColorMap.insert(pair<string, string>("Ac",        "70ABFA"));
    this->ColorMap.insert(pair<string, string>("Th",        "00BAFF"));
    this->ColorMap.insert(pair<string, string>("Pa",        "00A1FF"));
    this->ColorMap.insert(pair<string, string>("U",         "008FFF"));
    this->ColorMap.insert(pair<string, string>("Np",        "0080FF"));
    this->ColorMap.insert(pair<string, string>("Pu",        "006BFF"));
    this->ColorMap.insert(pair<string, string>("Am",        "545CF2"));
    this->ColorMap.insert(pair<string, string>("Cm",        "785CE3"));
    this->ColorMap.insert(pair<string, string>("Bk",        "8A4FE3"));
    this->ColorMap.insert(pair<string, string>("Cf",        "A136D4"));
    this->ColorMap.insert(pair<string, string>("Es",        "B31FD4"));
    this->ColorMap.insert(pair<string, string>("Fm",        "B31FBA"));
    this->ColorMap.insert(pair<string, string>("Md",        "B30DA6"));
    this->ColorMap.insert(pair<string, string>("No",        "BD0D87"));
    this->ColorMap.insert(pair<string, string>("Lr",        "C70066"));
    this->ColorMap.insert(pair<string, string>("Rf",        "CC0059"));
    this->ColorMap.insert(pair<string, string>("Db",        "D1004F"));
    this->ColorMap.insert(pair<string, string>("Sg",        "D90045"));
    this->ColorMap.insert(pair<string, string>("Bh",        "E00038"));
    this->ColorMap.insert(pair<string, string>("Hs",        "E6002E"));
    this->ColorMap.insert(pair<string, string>("Mt",        "EB0026"));
    this->ColorMap.insert(pair<string, string>("Ds",        "EB0026"));
    this->ColorMap.insert(pair<string, string>("Cn",        "EB0026"));
    this->ColorMap.insert(pair<string, string>("Nh",       "EB0026"));
    this->ColorMap.insert(pair<string, string>("Fl",        "EB0026"));
    this->ColorMap.insert(pair<string, string>("Mc",       "EB0026"));
    this->ColorMap.insert(pair<string, string>("Lv",        "EB0026"));
    this->ColorMap.insert(pair<string, string>("Og",       "EB0026"));
    this->ColorMap.insert(pair<string, string>("Ts",       "EB0026"));
    this->ColorMap.insert(pair<string, string>("00",        "98f58e"));

    // Constants and Conversions were taken from Gaussian'16 manual (https://gaussian.com/constants/)
    this->ConstantMap.insert(pair<string, double>("a0_bohr", 0.52917721092)); // in Angstrom
    this->ConstantMap.insert(pair<string, double>("qe_c", 1.602176565e-19)); // in Coulombs
    this->ConstantMap.insert(pair<string, double>("qe_esu", 4.803204e-10)); // in ESU
    this->ConstantMap.insert(pair<string, double>("h", 6.62606957e-34)); // in Joule-secs
    this->ConstantMap.insert(pair<string, double>("PI", 3.14159265358979323846)); // dimensionless
    this->ConstantMap.insert(pair<string, double>("hbar", 1.054571726e-34)); // in Joule-secs
    this->ConstantMap.insert(pair<string, double>("NA", 6.02214129e23)); // dimensionless
    this->ConstantMap.insert(pair<string, double>("hartree", 4.35974434e-18)); // in Joules
    this->ConstantMap.insert(pair<string, double>("C", 2.99792458e10)); // in cm/sec
    this->ConstantMap.insert(pair<string, double>("k_B", 1.3806488e-23)); // in Joule/Kelvin
    this->ConstantMap.insert(pair<string, double>("me", 9.10938291e-31)); // in kg
    this->ConstantMap.insert(pair<string, double>("R_cal", 1.9872036)); // in cal/(mol*K)
    this->ConstantMap.insert(pair<string, double>("R_J", 8.3144621)); // in J/(mol*K)
    this->ConstantMap.insert(pair<string, double>("P_0", 100000)); // in Pa

    this->ConversionMap.insert(pair<string, double>("Hartree_to_ev", 27.21138602)); // Hartree to eV conversion
    this->ConversionMap.insert(pair<string, double>("Hartree_to_KJ", 2625.499638)); // Hartree to KJ/mol conversion
    this->ConversionMap.insert(pair<string, double>("Hartree_to_Kcal", 627.509474)); // Hartree to Kcal/mol conversion
    this->ConversionMap.insert(pair<string, double>("ev_to_Joule", 1.602176565e-19)); // eV to Joule conversion
    this->ConversionMap.insert(pair<string, double>("Kcal_to_KJ", 4.184)); // Kcal to KJ conversion
    this->ConversionMap.insert(pair<string, double>("atm_to_Pa", 101325)); // atm to Pa conversion
    this->ConversionMap.insert(pair<string, double>("bar_to_Pa", 100000)); // bar to Pa conversion
    this->ConversionMap.insert(pair<string, double>("amu_to_kg", 1.660538921E-27)); // AMU to Kg conversion
    this->ConversionMap.insert(pair<string, double>("amu_to_kg_m2", 4.65082513926e-48)); // AMU to Kg/m^2 conversion

};

PeriodicTable::~PeriodicTable(){
    this->symbolMap.clear();
    this->massMap.clear();
    this->radiiMap.clear();
    this->ColorMap.clear();
    this->ConstantMap.clear();
    this->ConversionMap.clear();
};   

int PeriodicTable::getAtomicNumber(string symbol){
    return this->symbolMap[symbol];
};

double PeriodicTable::getConstant(string constantName){
    if (this->ConstantMap.find(constantName) == this->ConstantMap.end()){
        std::string msg = "Constant '" + constantName + "' not found. Available constants: ";
        bool first = true;
        for (const auto &p : this->ConstantMap){
            if (!first) msg += ", ";
            msg += p.first;
            first = false;
        }
        throw std::runtime_error(msg);
    };

    return this->ConstantMap[constantName];
};

double PeriodicTable::getConversion(string conversionName){
    if (this->ConversionMap.find(conversionName) == this->ConversionMap.end()){
        std::string msg = "Conversion '" + conversionName + "' not found. Available conversions: ";
        bool first = true;
        for (const auto &p : this->ConversionMap){
            if (!first) msg += ", ";
            msg += p.first;
            first = false;
        }
        throw std::runtime_error(msg);
    };
    return this->ConversionMap[conversionName];
};

double PeriodicTable::getAtomicMass(string symbol){
    return this->massMap[symbol];
};

double PeriodicTable::getCovalentRadii(string symbol){
    return this->radiiMap[symbol];
};

string PeriodicTable::getColor(string symbol){
    return this->ColorMap[symbol];
};

string PeriodicTable::getSymbol(int atomicNumber){
    string target = "";
    for (map<string, int>::iterator it=this->symbolMap.begin(); it!=this->symbolMap.end(); ++it){
        if (it->second == atomicNumber){
            return it->first;
        };
    };
    return target;

};

