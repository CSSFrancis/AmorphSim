from unittest import TestCase
import numpy as np
import matplotlib.pyplot as plt

from AmorphSim.utils.vector_utils import rotation_matrix_from_vectors, build_ico, build_ico_positions, get_ico_edges, get_ico_faces
from AmorphSim.clusters import Icosahedron
from mpl_toolkits.mplot3d import Axes3D

class TestVectorUtils(TestCase):

    def test_rot_matrix(self):
        rot = rotation_matrix_from_vectors([1,0,0],[1,1,1])
        print(np.dot(rot,[0,0,-4]))

    def test_build_icos(self):
        k,j,l = build_ico()
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(k[:,0],k[:,1],k[:,2], label="5-fold Axes")
        ax.scatter(j[:, 0], j[:, 1], j[:, 2],label="3-fold Axes")
        ax.scatter(l[:, 0], l[:, 1], l[:, 2],label="2-fold Axes")
        plt.legend()
        plt.show()

    def test_icoso_positions(self):
        pos = build_ico_positions()
        M = rotation_matrix_from_vectors(vec1=[0,0,1],vec2=pos[1])
        rot_pos = np.array([np.dot(p,M) for p in pos])
        for p in rot_pos:
            print(40, p[0]*2.7+100, p[1]*2.7+100, p[2]*2.7+20, 1.0, .8)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(rot_pos[:,0],rot_pos[:,1],rot_pos[:,2], s=500)
        plt.show()

    def test_get_ico_edges(self):
        i = Icosahedron("Zr", "Zr")
        atoms, edges = get_ico_edges(i.xyz[1:], shell=1)
        print(len(atoms))

    def test_get_ico_faces(self):
        i = Icosahedron("Zr", "Zr")
        faces, vectors = get_ico_faces(i.xyz[1:], shell=1)
        print(faces)
        print(vectors)