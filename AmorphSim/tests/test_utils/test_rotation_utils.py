from unittest import TestCase
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from AmorphSim.utils.rotation_utils import _get_rotation_matrix,  _get_random_rotation, _get_points_on_sphere
from AmorphSim.utils.simulation_utils import _shape_function, s_g_kernel



class TestRotUtils(TestCase):

    def test_get_random_rotation(self):
        rot, angles = _get_random_rotation(100)
        self.assertAlmostEqual(max(np.linalg.norm(rot,axis=1)), 1)
        self.assertAlmostEqual(min(np.linalg.norm(rot,axis=1)), 1)
        #fig = plt.figure()
        #ax = fig.add_subplot(111, projection='3d')
        #ax.scatter(rot[:,0],rot[:,1],rot[:,2])
        #plt.show()

    def test_random_rot_mat(self):
        rot = []
        for i in range(0,2000):
            M = rand_rotation_matrix(1.0)
            rot.append(M.dot([0,0,1]))
        rot= np.array(rot)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(rot[:,0],rot[:,1],rot[:,2])
        plt.show()

    def test_get_rotation_matrix(self):
        matrix = _get_rotation_matrix(axis=(1,1,0),theta=np.pi/3)
        print(matrix)
        print(np.dot(matrix,[4,0,0]))
        print(np.dot(matrix, [-4, 0, 0]))

    def test_get_deviation(self):
        sphere_rad = 1/_get_wavelength(200)
        matrix = _get_rotation_matrix(axis=(0,.2, .8), theta=np.pi / 3)
        k1=np.dot(matrix, [4, 0, 0])
        k2=np.dot(matrix, [-4, 0, 0])
        print(k1)
        dev1 =_get_deviation(sphere_radius=sphere_rad,k=k1)
        dev2 = _get_deviation(sphere_radius=sphere_rad,k=k2)
        print(dev1,dev2)
        int1=_shape_function(radius=.1, deviation=dev1)
        int2= _shape_function(radius=.1, deviation=dev2)
        print(int1,int2)

    def test_points_on_sphere(self):
        points = _get_points_on_sphere()
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(points[:,0],points[:,1],points[:,2])
        plt.show()