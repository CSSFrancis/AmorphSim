from unittest import TestCase
from AmorphSim.clusters import Icosahedron, FCC,  BCC
import matplotlib.pyplot as plt
import numpy as np


class TestIcosahedron(TestCase):
    def test_init(self):
        i = Icosahedron(central_atom="Zr", outer_atoms="Cu")
        i2 = Icosahedron(central_atom="Zr", outer_atoms="Cu", size=3)
        i3 = Icosahedron(central_atom="Zr", outer_atoms="Cu", size=3, anti=True)
        i3.plot()


class TestFCC(TestCase):
    def test_init(self):
        f = FCC(atom1="Cu", atom2="Cu", lattice_parameter=3.59, size=2.6)
        assert len(f) == 7

    def test_6fold_axes(self):
        f = FCC(atom1="Cu", atom2="Cu", lattice_parameter=3.59, size=5)
        f.get_6_fold_axis(inplace=True)
        f.plot()
        plt.show()


    def test_fcc(self):
        f = FCC(atom1="Zr", atom2="Cu", size=3)

    def test_bcc(self):
        b = BCC(atom1="Zr", atom2="Cu", size=3)
        b.plot(save=True)
        plt.show()

    def test_plot(self):
        b = BCC(atom1="Zr", atom2="Cu", size=1)
        b.plot(save=True)
        plt.show()

    def test_4foldFCC(self):
        b = FCC(atom1="Zr", atom2="Cu",size=1)
        #b.rotate_from_vectors(vector1=[0, 0, 1], vector2=[1, 1, 0], inplace=True)
        b.plot(save=True)
        plt.show()

    def test_6foldFCC(self):
        b = FCC(atom1="Zr", atom2="Cu", size=1)
        b.rotate_from_vectors(vector1=[0, 0, 1], vector2=[1, 1, 0], inplace=True)
        b.plot(save=True)
        plt.show()

    def test_4foldBCC(self):
        b = BCC(atom1="Zr", atom2="Cu",size=2)
        #b.rotate_from_vectors(vector1=[0, 0, 1], vector2=[0.0, -0.52573111, 0.85065081], inplace=True)
        b.plot(save=True)
        plt.show()

    def test_3foldBCC(self):
        b = BCC(atom1="Zr", atom2="Cu", size=2)
        b.rotate_from_vectors(vector1=[0, 0, 1], vector2=[1, 1, 1], inplace=True)
        b.plot(save=True)
        plt.show()


    def test_5fold(self):
        b = Icosahedron(central_atom="Zr", outer_atoms="Cu", size=1, anti=False)
        b.rotate_from_vectors(vector1=[0, 0, 1], vector2=[0.0, -0.52573111, 0.85065081], inplace=True)
        b.plot(save=True)
        plt.show()

    def test_3fold(self):
        b = Icosahedron(central_atom="Zr", outer_atoms="Cu", size=1, anti=False)
        b.rotate_from_vectors(vector1=[0, 0, 1], vector2=[-0.45879397, -0.45879397, 0.45879397], inplace=True)
        b.plot(save=True)
        plt.show()

    def test_2fold(self):
        b = Icosahedron(central_atom="Zr", outer_atoms="Cu", size=1, anti=False)
        #b.rotate_from_vectors(vector1=[0, 0, 1], vector2=[0.0, -0.52573111, 0.85065081], inplace=True)
        b.plot(save=True)
        plt.show()

    def test_all_rot(self):
        anti_mackay = AntiMackay(central_atom="Cu", outer_atoms="Zr")
        anti_mackay.all_direction_xyz(100)

