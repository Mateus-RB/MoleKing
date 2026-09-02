from MoleKing import G16LOGfile, ORCALOGfile, Psi4OUTfile
import math
import os
import platform

import pytest

class TestG16Output:
    @classmethod
    def setup_class(cls):
        cls.home_path = os.path.abspath(os.path.dirname(__file__))
        print(f"Home path: {cls.home_path}")

    def test_exist(self):
        assert os.path.exists(os.path.join(self.home_path, 'MK_test1.log')), "Log file does not exist"
        assert os.path.exists(os.path.join(self.home_path, 'MK_orbitals.log')), "Orbitals file does not exist"
        assert os.path.exists(os.path.join(self.home_path, 'MK_link.log')), "Link file does not exist"
        assert os.path.exists(os.path.join(self.home_path, 'MK_link2.log')), "Link2 file does not exist"
        assert os.path.exists(os.path.join(self.home_path, 'MK_polar.log')), "Polar file does not exist"

    def test_reader(self):
        log_path = os.path.join(self.home_path, 'MK_test1.log')
        file = G16LOGfile(log_path).__str__()
        if platform.system() == 'Windows':
            assert file == "G16LOGFile: Calculation of {} done in 2024, with the level of theory B3LYP/6-31+G(d) (6D, 7F) and SCF energy of -155.045622 Hartrees.".format(log_path)
        else:
            assert file == "G16LOGFile: Calculation of MK_test1.log done in 2024, with the level of theory B3LYP/6-31+G(d) (6D, 7F) and SCF energy of -155.045622 Hartrees."

    def test_getMol(self):
        Mol = G16LOGfile(os.path.join(self.home_path, 'MK_test1.log')).getMolecule()
        for atom in Mol:
            assert atom.getAtomicSymbol() in ['H', 'C', 'O']

    def test_getDipole(self):
        Ob = G16LOGfile(os.path.join(self.home_path,'MK_test1.log')).getDipole()
        assert Ob == 1.8039

    def test_getOrbitals(self):
        Ob = G16LOGfile(os.path.join(self.home_path,'MK_orbitals.log')).getOrbitals()
        assert len(Ob['Occupied']) == 13
        assert len(Ob['Unoccupied']) == 56

    def test_HOMO_LUMO(self):       
        Ob = G16LOGfile(os.path.join(self.home_path,'MK_orbitals.log'))
        assert Ob.getHOMO() == -0.27814
        assert Ob.getHOMO(-2) == -0.36919

        assert Ob.getLUMO() == 0.00034
        assert Ob.getLUMO(+4) == 0.07421

    def test_DetectLinks(self):
        Ob1 = G16LOGfile(os.path.join(self.home_path,'MK_link.log'), link=1)
        Ob2 = G16LOGfile(os.path.join(self.home_path,'MK_link.log'), link=2)

        assert Ob1.getEnergy() == -155.045622249
        assert Ob2.getEnergy() == -154.96646388

    def test_LinkOrbitals(self):
        Ob1 = G16LOGfile(os.path.join(self.home_path,'MK_link2.log'), link=1).getOrbitals()
        Ob2 = G16LOGfile(os.path.join(self.home_path,'MK_link2.log'), link=2).getOrbitals()      

        assert len(Ob1['Occupied']) == 13
        assert len(Ob2['Occupied']) == 5

    def test_getAlpha(self):
        f = G16LOGfile(os.path.join(self.home_path,'MK_polar.log'), polarAsw=True).getNLOFrequency()
        a1 = G16LOGfile(os.path.join(self.home_path,'MK_polar.log'), polarAsw=True).getAlpha(unit='esu', frequency=f[0])['iso']
        a2 = G16LOGfile(os.path.join(self.home_path,'MK_polar.log'), polarAsw=True).getAlpha(unit='SI', frequency=f[1])['xx']

        assert f == [0.0, 656.3, 587.6, 486.1]
        assert a1 == 95.8572
        assert a2 == 103.56

    def test_getBeta(self):
        f = G16LOGfile(os.path.join(self.home_path,'MK_polar.log'), polarAsw=True).getNLOFrequency()
        b1 = G16LOGfile(os.path.join(self.home_path,'MK_polar.log'), polarAsw=True).getBeta(unit='esu', frequency=f[0])['||']
        b2 = G16LOGfile(os.path.join(self.home_path,'MK_polar.log'), polarAsw=True).getBeta(unit='au', frequency=f[1], BSHG=1)['_|_(z)']

        assert b1 == 0.3999
        assert b2 == -0.0219886

    def test_getGamma(self):
        f = G16LOGfile(os.path.join(self.home_path,'MK_polar.log'), polarAsw=True).getNLOFrequency()
        g1 = G16LOGfile(os.path.join(self.home_path,'MK_polar.log'), polarAsw=True).getGamma(unit='esu', frequency=f[0])['||']
        g2 = G16LOGfile(os.path.join(self.home_path,'MK_polar.log'), polarAsw=True).getGamma(unit='au', frequency=f[1], GSHG=1)['xxxx']

        assert g1 == 31.0094
        assert g2 == 61491.6

    def test_thermo(self):
        f = G16LOGfile(os.path.join(self.home_path,'MK_Thermo.log'), thermoAsw=True)
        freq = f.getVibFrequencies()[-1]
        zpe = f.getZPE()
        zpve = f.getZPVE()
        h = f.getH()
        s = f.getS()
        g = f.getG()

        assert freq == 3864.4787
        assert zpe == -76.451023
        assert zpve == 13.20418
        assert h == -76.447243
        assert s == 45.112
        assert g == -76.468677

    def test_partition(self):
        f = G16LOGfile(os.path.join(self.home_path,'MK_Thermo.log'), thermoAsw=True)
        ele = abs((f.get_qEle() - 0.100000E+01)/ 0.100000E+01)
        vib = abs((f.get_qVib() - 0.100047E+01)/ 0.100047E+01)
        rot = abs((f.get_qRot() - 0.438927E+02)/ 0.438927E+02)
        trans = abs((f.get_qTrans() - 0.300431E+07)/ 0.300431E+07)
        
        assert ele < 0.05
        assert vib < 0.05
        assert rot < 0.05
        assert trans < 0.05


class TestORCAOutput:
    @staticmethod
    def write_thermo_log(
        path,
        atoms,
        rotational_constants,
        symmetry_number=1,
        frequency_blocks=((1000.0,),),
    ):
        """Write the minimum ORCA sections required by the output parser."""
        geometry = "\n".join(
            f"{symbol} {x:.6f} {y:.6f} {z:.6f}"
            for symbol, x, y, z in atoms
        )
        frequencies = "".join(
            "Scaling factor for frequencies =  1.000000000  (already applied!)\n"
            "----------------------------\n"
            + "\n".join(
                f"{index}: {frequency:.6f}"
                for index, frequency in enumerate(block)
            )
            + "\n\n"
            for block in frequency_blocks
        )
        bx, by, bz = rotational_constants
        path.write_text(
            f"""CARTESIAN COORDINATES (ANGSTROEM)
----------------------------
{geometry}
----------------------------
CARTESIAN COORDINATES (A.U.)
 Total Charge           Charge          ....    0
 Multiplicity           Mult            ....    1
 FINAL SINGLE POINT ENERGY                 -75.000000
{frequencies}
Point Group:  C1, Symmetry Number:   {symmetry_number}
Rotational constants in cm-1:     {bx:.6f}     {by:.6f}     {bz:.6f}
""",
            encoding="utf-8",
        )

    def test_linear_rotational_partition_function(self, tmp_path):
        log_path = tmp_path / "linear_orca.log"
        self.write_thermo_log(
            log_path,
            [("O", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 0.97)],
            (0.0, 18.774324, 18.774324),
        )

        output = ORCALOGfile(str(log_path), thermoAsw=True)
        theta_r = 18.774324 * 1.4387769599838156
        expected = 298.15 / theta_r

        # Linear formula from Thermochemistry in Gaussian, section 2.3:
        # q_rot = T / (sigma_r * Theta_r).
        assert math.isclose(output.get_qRot(), expected, rel_tol=1.0e-12)

    def test_nonlinear_rotational_partition_function(self, tmp_path):
        log_path = tmp_path / "nonlinear_orca.log"
        constants = (1.165006, 0.304730, 0.265422)
        self.write_thermo_log(
            log_path,
            [
                ("O", 0.0, 0.0, 0.0),
                ("H", 0.8, 0.0, 0.6),
                ("H", -0.8, 0.0, 0.6),
            ],
            constants,
            symmetry_number=2,
        )

        output = ORCALOGfile(str(log_path), thermoAsw=True)
        theta_x, theta_y, theta_z = [
            value * 1.4387769599838156 for value in constants
        ]
        expected = (
            math.sqrt(math.pi)
            * 298.15 ** 1.5
            / (2.0 * math.sqrt(theta_x * theta_y * theta_z))
        )

        # Nonlinear formula from the same manual section uses all three
        # rotational temperatures and the rotational symmetry number.
        assert math.isclose(output.get_qRot(), expected, rel_tol=1.0e-12)

    def test_qrot_requires_thermochemistry_data(self, tmp_path):
        log_path = tmp_path / "orca_without_thermo_flag.log"
        self.write_thermo_log(
            log_path,
            [("O", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 0.97)],
            (0.0, 18.774324, 18.774324),
        )

        output = ORCALOGfile(str(log_path))

        # Missing optional parsing must become a Python exception, never a
        # C++ out-of-bounds access/segmentation fault.
        with pytest.raises(RuntimeError, match="thermoAsw=True"):
            output.get_qRot()

    def test_uses_only_last_complete_frequency_block(self, tmp_path):
        log_path = tmp_path / "opt_ts_with_initial_and_final_hessians.log"
        self.write_thermo_log(
            log_path,
            [("O", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 0.97)],
            (0.0, 18.774324, 18.774324),
            frequency_blocks=(
                (-1690.70, 900.0),   # Initial Calc_Hess at the TS guess.
                (-1754.88, 1000.0),  # Final Freq at the optimized TS.
            ),
        )

        output = ORCALOGfile(str(log_path), thermoAsw=True)

        # Rate calculations must use the Hessian of the final optimized
        # geometry, never concatenate it with the initial OptTS Hessian.
        assert output.getVibFrequencies() == [-1754.88, 1000.0]


class TestPSI4Output():

    @classmethod
    def setup_class(cls):
        cls.home_path = os.path.abspath(os.path.dirname(__file__))

    def test_psi4(self): 
        file = Psi4OUTfile(os.path.join(self.home_path,'MK_Test.out'))
        assert file.getMul() == 1
        assert file.getCharge() == 0

    def test_psi4_geo(self):
        file = Psi4OUTfile(os.path.join(self.home_path,'MK_Test.out'))
        assert file.getMolecule().__str__() == 'Molecule BR_{4}O_{2}N_{2}C_{44}H_{30}, with charge 0 and multiplicity 1'
        assert file.getCharge() == 0
