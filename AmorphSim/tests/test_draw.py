from unittest import TestCase
from AmorphSim.clusters import FCC, Icosahedron
from AmorphSim.real_sim import Cube
from AmorphSim.draw.draw_3d import Icosohedron_3d, FCC_3d
import matplotlib.pyplot as plt
import numpy as np


class TestIco(TestCase):
    def test_draw_ico(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlim([-1, 1])
        ax.set_ylim([-1, 1])
        ax.set_zlim([-1, 1])
        ax.axis("off")
        i = Icosohedron_3d(line_width=1, alpha=.5,radius=1)
        ax.add_collection3d(i)
        plt.show()

    def test_draw_fcc(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlim([-3, 3])
        ax.set_ylim([-3, 3])
        ax.set_zlim([-3, 3])
        i = FCC_3d(line_width=1, alpha=.5)
        ax.add_collection3d(i)
        from itertools import permutations
        a = [0, 1, 2]
        perms = list(permutations(a))
        a = [0, -1, 2]
        perms = perms +list(permutations(a))
        a = [0, 1, -2]
        perms = perms+(list(permutations(a)))
        a = [0, -1, -2]
        perms = perms + list(permutations(a))
        perms = np.array(perms)
        print(perms)
        ax.scatter(perms[:,0], perms[:,1], perms[:,2])
        plt.show()

    def test_draw_cube(self):
        i = Icosahedron(position=[14, 3, 10], num_shells=1)
        i2 = Icosahedron(position=[10, 10, 3], num_shells=1)
        i3 = Icosahedron(position=[3, 17, 16], num_shells=1)
        i4 = Icosahedron(position=[13, 17, 16], num_shells=1)
        i5 = Icosahedron(position=[3, 3, 5], num_shells=1)
        f = FCC(position=[1, 13, 6], radius=2)
        f2 = FCC(position=[12, 3, 2], radius=1.5)
        f3 = FCC(position=[12, 17, 2], radius=1.5)
        f4 = FCC(position=[12, 13, 16], radius=1.5)
        f5 = FCC(position=[11, 14, 5], radius=1.5)
        c = Cube()
        c.add_cluster(i)
        c.add_cluster(i2)
        c.add_cluster(i3)
        c.add_cluster(i4)
        c.add_cluster(i5)
        c.add_cluster(f)
        c.add_cluster(f2)
        c.add_cluster(f3)
        c.add_cluster(f4)
        c.add_cluster(f5)
        c.plot_3d()
        plt.show()
