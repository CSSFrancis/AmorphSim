from unittest import TestCase
import numpy as np
import matplotlib.pyplot as plt
import hyperspy.api as hs
from AmorphSim.sim import Cluster, SimulationCube
from AmorphSim.utils.rotation_utils import _rand_2d_rotation_matrix
from mpl_toolkits.mplot3d import Axes3D


class TestCluster(TestCase):
    def setUp(self):
        self.c = Cluster(symmetry=10,radius=.5, k=4.0, position=(1,1))

    def test_get_k_vectors(self):

        self.c.rotation_2d = _rand_2d_rotation_matrix()
        print(self.c.rotation_2d)
        k = np.array(self.c.get_k_vectors())
        plt.scatter(k[:,0], k[:,1])
        plt.show()
    def test_get_k_vectors_rot(self):
        #self.c.rotation_3d = np.eye(3)
        self.c.plane_direction=[1,1,1]
        print(self.c.rotation_2d)
        k = np.array(self.c.get_k_vectors())
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(k[:,0],k[:,1],k[:,2])
        self.c.plane_direction = [1, -1, 1]
        print(self.c.rotation_2d)
        k = np.array(self.c.get_k_vectors())
        ax.scatter(k[:, 0], k[:, 1], k[:, 2])
        plt.show()

    def test_get_diffraciton(self):
        diff = self.c.get_diffraction(img_size=15.0)
        plt.imshow(diff)
        plt.show()

class TestSimulationCube(TestCase):
    def test_random_init(self):
        cube = SimulationCube()
        cube.add_random_clusters(100)
        print(cube)

    def test_ico_init(self):
        cube = SimulationCube()
        cube.add_icoso(1, radius_range=(4., 4.1))
        stem = cube.get_4d_stem(noise=True,convolve=True)
        stem.plot()
        plt.show()
        print(cube)

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



