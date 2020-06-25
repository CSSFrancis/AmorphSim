from unittest import TestCase
from AmorphSim.real_sim import Cluster
from AmorphSim.clusters import Cubic,Icosohedron,AntiMckay,FCC,BCC


class TestRealCube(TestCase):
    def test_to_xyz(self):
        b = BCC(atom1="Zr", atom2="Cu")
        b.get_xyz()