from unittest import TestCase
from AmorphSim.clusters import BCC


class TestRealCube(TestCase):
    def test_to_xyz(self):
        b = BCC(atom1="Zr", atom2="Cu")
        b.get_xyz()