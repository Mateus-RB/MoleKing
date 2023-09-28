from MoleKing import Molecule

class TestMolecule():

    def test_add(self):
        mol = Molecule()
        mol.addAtom('H', 0, 0, 0)
        mol.addAtom('H', 0, 0, 1)
        for atom in mol:
            assert atom.getAtomicSymbol() == 'H'

    def test_len(self):
        mol = Molecule()
        mol.addAtom('H', 0, 0, 0)
        mol.addAtom('H', 0, 0, 1)
        assert len(mol) == 2 

    def test_getBonds(self):
        mol = Molecule()
        mol.addAtom('O', 0, 0, 0)
        mol.addAtom('H', 0, 0, 1)
        mol.addAtom('H', 0, 1, 0)
        assert [0, 1] in mol.getIRCBonds()

    def test_getangles(self):
        mol = Molecule()
        mol.addAtom('O', 0, 0, 0)
        mol.addAtom('H', 0, 0, 1)
        mol.addAtom('H', 0, 1, 0)
        assert [1, 0, 2] in mol.getIRCAngles()

    def test_getdihedrals(self):
        mol = Molecule()
        mol.addAtom('C', 0, 0, 0)
        mol.addAtom('H', 0, 0, 1)
        mol.addAtom('H', 0, 1, 0)
        mol.addAtom('O', 1, 0, 0)
        mol.addAtom('H', 2, 0, 0)
        assert [1, 0, 3, 4] in mol.getIRCDihedrals()

    def test_center(self):
        mol = Molecule()
        mol.addAtom('O', 0, 0, 0)
        mol.addAtom('H', 0, 0, 1)
        mol.addAtom('H', 0, 1, 0)
        C = mol.getMassCenter()
        assert C.getCoords('c')[1] < 0.06 
    
    def test_translate(self):
        mol = Molecule()
        mol.addAtom('O', 0, 0, 0)
        mol.addAtom('H', 0, 0, 1)
        mol.addAtom('H', 0, 1, 0)
        mol.moveMassCenter(10, 10, 10)
        C = mol.getMassCenter()
        assert C.getCoords('c') == [10.0, 10.0, 10.0]
    
    def test_VDW(self):
        mol = Molecule()
        mol.setVDWRatio(2)
        assert mol.getVDWRatio() == 2
