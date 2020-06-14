from builtins import len

import numpy as np
from diffpy.structure import Atom, Structure
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

class Cube:
    """
    The main idea with the Cube object is that it is a simplified way of visualizing a glass as a list
    of clusters and their relative positions. From that idea we can build a model of each cluster
    as well as get the positions of each atom for simulation either dynamically or not...
    """
    def __init__(self, dimensions=(20, 20, 20)):
        """Initializes the simulation cube in nm
        Parameters
        --------------
        dimensions: tuple
            The dimensions of the cube to simulate.
        """
        self.clusters = []
        self.dimensions = np.array(dimensions)

    def __str__(self):
        return ("<Cube of " + len(self.clusters) + " clusters")

    def to_prismatic_xyz(self):
        """Prints to format for simulation with Prismatic. (.xyz files)
        """
        pass

    def prism_simulate(self,beam_direction=[1,0,0], **kwargs):
        """Simulates using Prismatic to return a 4-D STEM dataset along some beam direction...
        """
        pass
    def kinematic_simulate(self,
                           beam_size=1,
                           beam_step=1.,
                           beam_direction=[1,0,0],
                           ):
        """Simulates using kinematic diffraction
        :param beam_direction:
        :param kwargs:
        :return:
        """
        pass

    def plot_2d(self):
        """Plots the clusters in 2 dimensional projection
        """
        pass

    def plot_3d(self):
        """Plots the clusters in 3 dimensions as prototypical shapes?
        """
        pass


class Cluster(Structure):
    """Each Cluster extends the Structure class giving some unique simulation abilities to the strucutre
    """
    def plot(self):
        """Plots the atoms of some structure in 3-D.  The atoms are shown as spheres based on their atomic
        size.
        """
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        coordinates = self.xyz_cartn
        print(self.element)
        ax.scatter(coordinates[:,0], coordinates[:,1], coordinates[:,2],
                   marker='o', s=300)


    def plot_reciprocal(self, resolution=100):
        """Plots the reciprocal space  atoms of some structure in 3-D. Gives a voxel representation of
        the space.
        """
        pass
    def get_reciprocal_space(self, resolution=100):
        """Returns the three dimensional reciprocal space fro some cluster where some cut along the
        reciprocal space represents the kinetic diffraction along some direction. Creates a resolution^3
        space
        """
        pass
