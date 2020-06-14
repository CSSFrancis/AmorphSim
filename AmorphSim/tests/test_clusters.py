from unittest import TestCase
from AmorphSim.clusters import Cubic,Icosohedron,Super_Icosohedron,FCC,BCC
import matplotlib.pyplot as plt
class TestClusters(TestCase):
    def test_icosahedron(self):
        i = Icosohedron(central_atom="Zr", outer_atoms="Cu")
        print(i.xyz_cartn[:,2])
        print()

    def test_plot(self):
        i = Icosohedron(central_atom="Zr", outer_atoms="Cu")
        i.plot()
        plt.show()
