from builtins import len

import numpy as np
import copy as copymod
from diffpy.structure import Atom, Structure,Lattice
from mpl_toolkits.mplot3d import Axes3D
from AmorphSim.utils.rotation_utils import _rand_2d_rotation_matrix,_rand_3d_rotation_matrix, _get_points_on_sphere
import os
from AmorphSim.utils.vector_utils import rotation_matrix_from_vectors,_get_angle_between
import matplotlib.pyplot as plt



class GlassCube:
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
        return "<Cube of " + len(self.clusters) + " clusters"

    def add_cluster(self, cluster):
        self.clusters.append(cluster)

    def to_prismatic_xyz(self):
        """Prints to format for simulation with Prismatic. (.xyz files)
        """
        pass

    def prism_simulate(self, beam_direction=[1, 0, 0], **kwargs):
        """Simulates using Prismatic to return a 4-D STEM dataset along some beam direction...
        """
        pass

    def kinematic_simulate(self,
                           beam_size=1,
                           beam_step=1.,
                           beam_direction=[1, 0, 0],
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
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlim([0, self.dimensions[0]])
        ax.set_ylim([0, self.dimensions[1]])
        ax.set_zlim([0, self.dimensions[2]])
        ax.set_xlabel("Real Space, X, nm")
        ax.set_ylabel("Real Space, Y, nm")
        #ax.set_xticks([])
        #ax.set_yticks([])
        ax.set_zticks([])
        ax.set_facecolor("grey")
        for c in self.clusters:
            d = c.draw()
            ax.add_collection3d(d)

