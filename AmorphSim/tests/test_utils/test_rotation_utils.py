from unittest import TestCase
import numpy as np

from AmorphSim.utils.rotation_utils import _get_rotation_matrix, _get_deviation, _get_wavelength,_shape_function


class TestRotUtils(TestCase):

    def test_get_rotation_matrix(self):
        matrix = _get_rotation_matrix(axis=(1,1,0),theta=np.pi/3)
        print(matrix)
        print(np.dot(matrix,[4,0,0]))
        print(np.dot(matrix, [-4, 0, 0]))

    def test_get_wavelength(self):
        self.assertAlmostEqual(_get_wavelength(200), 0.002508,4)
        self.assertAlmostEqual(_get_wavelength(300), 0.0019687, 4)

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