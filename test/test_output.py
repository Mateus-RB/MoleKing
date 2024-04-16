from MoleKing import G16LOGfile, Psi4OUTfile
from os import getcwd

class TestG16Output():

    def test_reader(self):
        home_path = getcwd()+'/test'
        if '/test/test' in home_path:
            home_path = home_path.replace('/test/test', '/test')    
        file = G16LOGfile(home_path+'/MK_test1.log').__str__()
        assert file == "G16LOGFile: Calculation of MK_test1.log done in 2024, with the level of theory B3LYP/6-31+G(d) (6D, 7F) and SCF energy of -155.045622 Hartrees."

    def test_getMol(self):
        home_path = getcwd()+'/test'
        if '/test/test' in home_path:
            home_path = home_path.replace('/test/test', '/test')    
        Mol = G16LOGfile(home_path+'/MK_test1.log').getMolecule()
        for atom in Mol:
            assert atom.getAtomicSymbol() in ['H', 'C', 'O']

    def test_getDipole(self):
        home_path = getcwd()+'/test'
        if '/test/test' in home_path:
            home_path = home_path.replace('/test/test', '/test')    
        Ob = G16LOGfile(home_path+'/MK_test1.log').getDipole()
        assert Ob == 1.8039

    def test_getOrbitals(self):
        home_path = getcwd()+'/test'
        if '/test/test' in home_path:
            home_path = home_path.replace('/test/test', '/test')    
        Ob = G16LOGfile(home_path+'/MK_orbitals.log').getOrbitals()
        assert len(Ob['Occupied']) == 13
        assert len(Ob['Unoccupied']) == 56

    def test_HOMO_LUMO(self):
        home_path = getcwd()+'/test'
        if '/test/test' in home_path:
            home_path = home_path.replace('/test/test', '/test')          
        Ob = G16LOGfile(home_path+'/MK_orbitals.log')
        assert Ob.getHOMO() == -0.27814
        assert Ob.getHOMO(-2) == -0.36919

        assert Ob.getLUMO() == 0.00034
        assert Ob.getLUMO(+4) == 0.07421

    def test_DetectLinks(self):
        home_path = getcwd()+'/test'
        if '/test/test' in home_path:
            home_path = home_path.replace('/test/test', '/test')    
        Ob1 = G16LOGfile(home_path+'/MK_link.log', link=1)
        Ob2 = G16LOGfile(home_path+'/MK_link.log', link=2)

        assert Ob1.getEnergy() == -155.045622249
        assert Ob2.getEnergy() == -154.96646388

    def test_LinkOrbitals(self):
        home_path = getcwd()+'/test'
        if '/test/test' in home_path:
            home_path = home_path.replace('/test/test', '/test')    
        Ob1 = G16LOGfile(home_path+'/MK_link2.log', link=1).getOrbitals()
        Ob2 = G16LOGfile(home_path+'/MK_link2.log', link=2).getOrbitals()      

        assert len(Ob1['Occupied']) == 13
        assert len(Ob2['Occupied']) == 5

    def test_getAlpha(self):
        home_path = getcwd()+'/test'
        if '/test/test' in home_path:
            home_path = home_path.replace('/test/test', '/test')   
        print(home_path+'/MK_polar.log')
        f = G16LOGfile(home_path+'/MK_polar.log', polarAsw=True).getFrequency()
        a1 = G16LOGfile(home_path+'/MK_polar.log', polarAsw=True).getAlpha(unit='esu', frequency=f[0])['iso']
        a2 = G16LOGfile(home_path+'/MK_polar.log', polarAsw=True).getAlpha(unit='SI', frequency=f[1])['xx']

        assert f == [0.0, 656.3, 587.6, 486.1]
        assert a1 == 95.8572
        assert a2 == 103.56

    def test_getBeta(self):
        home_path = getcwd()+'/test'
        if '/test/test' in home_path:
            home_path = home_path.replace('/test/test', '/test')   
        f = G16LOGfile(home_path+'/MK_polar.log', polarAsw=True).getFrequency()
        b1 = G16LOGfile(home_path+'/MK_polar.log', polarAsw=True).getBeta(unit='esu', frequency=f[0])['||']
        b2 = G16LOGfile(home_path+'/MK_polar.log', polarAsw=True).getBeta(unit='au', frequency=f[1], BSHG=1)['_|_(z)']

        assert b1 == 0.3999
        assert b2 == -0.0219886

    def test_getGamma(self):
        home_path = getcwd()+'/test'
        if '/test/test' in home_path:
            home_path = home_path.replace('/test/test', '/test')   
        f = G16LOGfile(home_path+'/MK_polar.log', polarAsw=True).getFrequency()
        g1 = G16LOGfile(home_path+'/MK_polar.log', polarAsw=True).getGamma(unit='esu', frequency=f[0])['||']
        g2 = G16LOGfile(home_path+'/MK_polar.log', polarAsw=True).getGamma(unit='au', frequency=f[1], GSHG=1)['xxxx']

        assert g1 == 31.0094
        assert g2 == 61491.6

class TestPSI4Output():

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
