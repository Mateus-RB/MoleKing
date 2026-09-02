from MoleKing import Molecule
import pytest

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

    @pytest.mark.parametrize(
        ("hf_percent", "expected_iop"),
        [
            (0, "IOP(3/76=1000000000)"),
            (5, "IOP(3/76=0950000500)"),
            (20, "IOP(3/76=0800002000)"),
            (100, "IOP(3/76=0000010000)"),
        ],
    )
    def test_toGJF_modHF(self, tmp_path, hf_percent, expected_iop):
        mol = Molecule()
        mol.addAtom("H", 0, 0, 0)
        output = tmp_path / f"modhf_{hf_percent}.gjf"

        mol.toGJF(fileName=str(output), method="M06HF", modHF=hf_percent)

        gjf = output.read_text()
        assert expected_iop in gjf
        assert "IOP(3/77=" not in gjf

    @pytest.mark.parametrize("method", ["B3LYP", "b3lyp", "RB3LYP", "UB3LYP"])
    def test_toGJF_modHF_preserves_b3lyp_exchange_ratio(
        self, tmp_path, method
    ):
        mol = Molecule()
        mol.addAtom("H", 0, 0, 0)
        output = tmp_path / f"{method}_modhf.gjf"

        mol.toGJF(
            fileName=str(output),
            method=method,
            modHF=20,
        )

        gjf = output.read_text()
        assert "IOP(3/76=0800002000)" in gjf
        assert "IOP(3/77=0900010000)" in gjf

    def test_toGJF_default_b3lyp_uses_b3lyp_exchange_ratio(self, tmp_path):
        mol = Molecule()
        mol.addAtom("H", 0, 0, 0)
        output = tmp_path / "default_b3lyp_modhf.gjf"

        mol.toGJF(fileName=str(output), modHF=20)

        gjf = output.read_text()
        assert "IOP(3/76=0800002000)" in gjf
        assert "IOP(3/77=0900010000)" in gjf

    @pytest.mark.parametrize("hf_percent", [-2, 101])
    def test_toGJF_rejects_invalid_modHF(self, tmp_path, hf_percent):
        mol = Molecule()
        mol.addAtom("H", 0, 0, 0)

        with pytest.raises(ValueError, match="modHF must be between 0 and 100"):
            mol.toGJF(
                fileName=str(tmp_path / f"invalid_{hf_percent}.gjf"),
                modHF=hf_percent,
            )
