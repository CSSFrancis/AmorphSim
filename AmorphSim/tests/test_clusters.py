from unittest import TestCase
from AmorphSim.clusters import Cubic, Icosohedron, Mackay, FCC,  BCC,AntiMackay
import matplotlib.pyplot as plt


class TestClusters(TestCase):
    def test_icosahedron(self):
        i = Icosohedron(central_atom="Zr", outer_atoms="Cu")
        print(i.xyz_cartn[:, 2])
        print()

    def test_bcc(self):
        b = BCC(atom1="Zr", atom2="Cu", size=4)
        print(len(b))
        b.get_xyz()

    def test_plot(self):
        b = BCC(atom1="Zr", atom2="Cu",size=4)
        b.plot()
        plt.show()

    def test_rot(self):
        b = BCC(atom1="Zr", atom2="Cu")
        b.rotate_from_vectors([1, 0, 0], [1, 1, 1])
        b.get_xyz(file="111.xyz")

    def test_random(self):
        b = Icosohedron(central_atom="Zr", outer_atoms="Cu")
        b.random_rot_xyz()

    def test_Mackay(self):
        mackay = Mackay(central_atom="Cu", outer_atoms="Zr")
        mackay.plot()
        plt.show()

    def test_anti_Mackay(self):
        mackay = AntiMackay(central_atom="Cu", outer_atoms="Zr")
        mackay.plot()
        plt.show()

