from MoleKing import G16LOGfile, Psi4OUTfile
from os import getcwd

class TestOutput():

    def test_reader(self):
        home_path = getcwd()+'/test'
        if '/test/test' in home_path:
            home_path = home_path.replace('/test/test', '/test')    
        file = G16LOGfile(home_path+'/MK_Test.log').__str__()
        assert file == "G16LOGFile: Calculation of MK_Test.log done in 16 at Mon Sep 25 14:59:53 2023, with the level of theory B3LYP/6-31+(d') (6D, 7F) and SCF energy of -155.043566 Hartrees."

    def test_getMol(self):
        home_path = getcwd()+'/test'
        if '/test/test' in home_path:
            home_path = home_path.replace('/test/test', '/test')    
        Mol = G16LOGfile(home_path+'/MK_Test.log').getMolecule()
        for atom in Mol:
            assert atom.getAtomicSymbol() in ['H', 'C', 'O']

    def test_getOrbitals(self):
        home_path = getcwd()+'/test'
        if '/test/test' in home_path:
            home_path = home_path.replace('/test/test', '/test')    
        Ob = G16LOGfile(home_path+'/MK_Test.log').getOrbitals()
        assert len(Ob['Occupied']) == 39
        assert len(Ob['Unoccupied']) == 168

    def test_getDipole(self):
        home_path = getcwd()+'/test'
        if '/test/test' in home_path:
            home_path = home_path.replace('/test/test', '/test')    
        Ob = G16LOGfile(home_path+'/MK_Test.log').getDipole()
        assert Ob == 1.8442

    def test_psi4(self):
        home_path = getcwd()+'/test'
        if '/test/test' in home_path:
            home_path = home_path.replace('/test/test', '/test')    
        file = Psi4OUTfile(home_path+'/MK_Test.out')
        assert file.getMul() == 1
        assert file.getCharge() == 0

    def test_psi4_geo(self):
        home_path = getcwd()+'/test'
        if '/test/test' in home_path:
            home_path = home_path.replace('/test/test', '/test')    
        file = Psi4OUTfile(home_path+'/MK_Test.out')
        assert file.getMolecule().__str__() == 'Molecule BR_{4}O_{2}N_{2}C_{44}H_{30}, with charge 0 and multiplicity 1'
        assert file.getCharge() == 0
