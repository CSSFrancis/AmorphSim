from unittest import TestCase
import numpy as np

from AmorphSim.utils.simulation_utils import _get_wavelength, _get_speckle_size, _get_deviation
from AmorphSim.utils.simulation_utils import _get_disorder

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