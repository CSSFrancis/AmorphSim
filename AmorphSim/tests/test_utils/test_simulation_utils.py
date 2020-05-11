from unittest import TestCase
import numpy as np
import matplotlib.pyplot as plt

from AmorphSim.utils.simulation_utils import _get_wavelength, _get_speckle_size, _get_deviation, _shape_function
from AmorphSim.utils.simulation_utils import _get_disorder, _get_speckle_intensity

class TestSimulationUtils(TestCase):

    def test_get_wavelength(self):
        self.assertAlmostEqual(_get_wavelength(200), 0.002508,4)
        self.assertAlmostEqual(_get_wavelength(300), 0.0019687, 4)

    def test_get_speckle_size(self):
        size = _get_speckle_size(200, .6)
        self.assertAlmostEqual(size, 0.23922, 4)

    def test_get_deviation(self):
        deviation = _get_deviation(100.0, [4,4,0])
        self.assertAlmostEqual(deviation,-0.15987220439131988)

    def test_get_disorder(self):
        f = _get_disorder([4,4,0],.01)
        print(f)

    def test_get_speckle_intensity(self):
        f = _get_speckle_intensity(k_vector=[4,4,0], ewald_sphere_rad=100000)
        print(f)

    def test_shape_function(self):
        a1,i1 = self.shape_function_vs_angle(1.0)
        a2, i2 = self.shape_function_vs_angle(2.0)
        a4, i4 = self.shape_function_vs_angle(4.0)
        plt.plot(np.degrees(a1),i1, label="1 nm cluster")
        plt.plot(np.degrees(a2), i2, label="2 nm cluster")
        plt.plot(np.degrees(a4), i4, label="4 nm cluster")
        plt.legend()
        plt.xlabel("Angle of Mistilt Degrees")
        plt.show()

    def shape_function_vs_angle(self, size=1.0):
        angles = np.linspace(start=-np.pi / 8, stop=np.pi / 8, num=100)
        rad = 1 / _get_wavelength(200)
        inten = []
        for angle in angles:
            M = np.array([[np.cos(angle), -np.sin(angle), 0],
                            [0, 0, 1],
                            [np.sin(angle), np.cos(angle), 0]
                            ])
            k = np.dot(M, [4, 0, 0])
            k[2] = -rad + k[2]
            disp = np.linalg.norm(k) - rad

            inten.append(_shape_function(radius=size, deviation=disp))
        return angles,inten