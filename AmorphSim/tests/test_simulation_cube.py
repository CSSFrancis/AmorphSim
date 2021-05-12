import pytest
import numpy as np
import matplotlib.pyplot as plt
import hyperspy.api as hs

import diffsims
from diffsims.generators.diffraction_generator import AtomicDiffractionGenerator
from AmorphSim.sim.simulation_cube import Cube
from AmorphSim.clusters.clusters import FCC,Icosahedron

class TestSimulationCube:

    @pytest.fixture
    def sim_cube(self):
        cube = Cube(dimensions=(200, 200, 20))
        cube.add_cluster(Icosahedron(position=(100, 100, 10)))
        cube.add_cluster(FCC(position=(120, 120, 10)))
        return cube

    def test_plot_3d(self, sim_cube):
        sim_cube.plot_3d()
        plt.show()

    def test_to_xyz(self, sim_cube):
        test_str = sim_cube.to_prismatic_xyz()
        print(test_str)
        test_str = sim_cube.to_prismatic_xyz(file="test_cube")

    def test_projection(self):
        cube = SimulationCube()
        cube.add_random_clusters(100)
        projection = cube.show_projection()
        plt.imshow(projection)
        plt.show()

    def test_symmetry_projection(self):
        cube = SimulationCube()
        cube.add_random_clusters(1000)
        projection = cube.plot_symmetries()

    def test_symmetry_projection_acceptance(self):
        cube = SimulationCube()
        cube.add_random_clusters(1000)
        projection = cube.plot_symmetries(acceptance=np.pi/8)
        #projection = cube.plot_symmetries(acceptance=np.pi / 2)

    def test_4dStem(self):
        cube = SimulationCube()
        cube.add_random_clusters(10000)
        #cube.plot_symmetries()
        stem = cube.get_4d_stem(noise=True, disorder=.05, convolve=True)
        stem.plot()
        plt.show()




"""def get_mistilt(n=10, radius=2, accept=np.pi/8,):
    theta = np.linspace(0, np.pi / 4, num=30)  # 22.5 deg
    phi = np.linspace(0, np.pi, num=30)[np.newaxis].T  #
    x, y, z = np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta) * np.ones((15, 1))
    power_list = []
    for phix, phiy, phiz in zip(x, y, z):
        pow_list = []
        for thetax, thetay, thetaz in zip(phix, phiy, phiz):
            c = Cluster(symmetry=n, radius=radius, k=4.0, position=(1, 1), rotation_vector=[thetax, thetay, thetaz],
                        rotation_angle=np.pi)
            print(thetax, thetay, thetaz)
            diff_image = c.get_diffraction(img_size=15.0)
            diff_image = diff_image + np.random.random(size=(512, 512))*10
            ps = np.sum(power_spectrum(angular_correlation(to_polar_image(diff_image, radius=[0, 200],
                                                                          center=(256, 256))))[:, n])
            pow_list.append(ps)
        power_list.append(pow_list)
    plt.errorbar(theta / np.pi * 180, np.mean(power_list, axis=0), yerr=np.std(power_list, axis=0),
                 label=str(n) + "-fold Symmetry")"""



