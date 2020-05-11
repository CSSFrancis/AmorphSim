from unittest import TestCase
import numpy as np
import matplotlib.pyplot as plt

from AmorphSim.utils.vector_utils import rotation_matrix_from_vectors,build_ico
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