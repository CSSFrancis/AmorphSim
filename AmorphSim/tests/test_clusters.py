from unittest import TestCase
from AmorphSim.clusters import Icosahedron, FCC,  BCC
import matplotlib.pyplot as plt
import numpy as np


class TestIcosahedron(TestCase):
    def test_init(self):
        i = Icosahedron(central_atom="Cu", outer_atoms="Cu", shell_distance=2.54, num_shells=1)
        assert isinstance(i, Icosahedron)
        assert len(i) == 13

    def test_mackay(self):
        i = Icosahedron(central_atom="Cu", outer_atoms="Cu", shell_distance=2.54, num_shells=2)
        assert isinstance(i, Icosahedron)
        assert len(i) == 55

    def test_antimackay(self):
        i = Icosahedron(central_atom="Cu", outer_atoms="Cu", shell_distance=2.54, num_shells=2, anti=True)
        assert isinstance(i, Icosahedron)
        assert len(i) == 45

    def test_get_2fold_axis(self):
        i = Icosahedron(central_atom="Cu", outer_atoms="Cu", shell_distance=2.54, num_shells=1)
        i.plot(rotate=False)
        plt.show()

    def test_get_5fold_axis(self):
        i = Icosahedron(central_atom="Cu", outer_atoms="Cu", shell_distance=2.54, num_shells=1)
        i.get_5_fold_axis(reinitialize=True)
        i.plot(rotate=False)
        plt.show()

    def test_get_3fold_axis(self):
        i = Icosahedron(central_atom="Cu", outer_atoms="Cu", shell_distance=2.54, num_shells=1)
        i.get_3_fold_axis(reinitialize=True)
        i.plot(rotate=False)
        plt.show()


class TestFCC(TestCase):
    def test_init(self):
        f = FCC(atom1="Cu", atom2="Cu", lattice_parameter=3.59, radius=2.6)
        assert len(f) == 7

    def test_6fold_axes(self):
        f = FCC(atom1="Cu", atom2="Cu", lattice_parameter=3.59, radius=5)
        f.get_6_fold_axis(inplace=True)
        f.plot(rotate=False)
        plt.show()

    def test_4fold_axes(self):
        f = FCC(atom1="Cu", atom2="Cu", lattice_parameter=3.59, radius=5)
        f.get_4_fold_axis(inplace=True)
        f.plot(rotate=False)
        plt.show()


class TestBCC(TestCase):
    def test_init(self):
        b = BCC(atom1="Fe", atom2="Fe", lattice_parameter=2.856, radius=3)
        assert len(b) == 15

    def test_6fold_axes(self):
        b = BCC(atom1="Fe", atom2="Fe", lattice_parameter=2.856, radius=5)
        b.get_6_fold_axis(inplace=True)
        b.plot(rotate=False)
        plt.show()

    def test_4fold_axes(self):
        f = BCC(atom1="Cu", atom2="Cu", lattice_parameter=2.856, radius=5)
        f.get_4_fold_axis(inplace=True)
        f.plot(rotate=False)
        plt.show()

class TestCluster(TestCase):
    def test_add_disorder(self):
        b = BCC(atom1="Cu", atom2="Cu", lattice_parameter=3.59, radius=5)
        b.add_disorder(sigma=.2)
        b.plot()
        plt.show()
