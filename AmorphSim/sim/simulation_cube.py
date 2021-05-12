from builtins import len

import numpy as np
import matplotlib.pyplot as plt
from diffpy.structure import Atom, Structure,Lattice
from AmorphSim.draw.draw_3d import Cube3d


class Cube(Structure):
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
        self.lattice = Lattice(a=1, b=1, c=1, alpha=90, beta=90, gamma=90)
        self.dimensions = np.array(dimensions)

    def __str__(self):
        return "<Cube of " + len(self.clusters) + " clusters"

    def add_cluster(self, cluster):
        """This function adds the atoms from some cluster to the
        :param cluster:
        :return:
        """
        for i in cluster:
            self.append(i)
        self.clusters.append(cluster)

    def to_prismatic_xyz(self,file=None,):
        """Prints to format for simulation with Prismatic. (.xyz files)
        """
        newstr = ("Simulating XYZ \n     " + str(self.dimensions[0]) + " " +
                  str(self.dimensions[1]) + " " + str(self.dimensions[2]) + " \n")

        for cluster in self.clusters:
            newstr = newstr + cluster.get_simple_xyz()

        if file is None:
            return newstr
        else:
            with open(file + ".xyz", "+w") as f:
                f.write(newstr)


    def prism_simulate(self, beam_direction=[1, 0, 0], **kwargs):
        """Simulates using Prismatic to return a 4-D STEM dataset along some beam direction...
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
        ax.set_xlim([-5, self.dimensions[0]+5])
        ax.set_ylim([-5, self.dimensions[1]+5])
        ax.set_zlim([-5, self.dimensions[0]+5])
        ax.set_xlabel("Real Space, X, nm")
        ax.set_ylabel("Real Space, Y, nm")
        #ax.set_xticks([])
        #ax.set_yticks([])
        ax.set_zticks([])
        #ax.set_facecolor("grey")
        c = Cube3d(shape=self.dimensions, alpha=0.2,)
        ax.add_collection3d(c)
        for c in self.clusters:
            d = c.draw()
            ax.add_collection3d(d)
        ax.grid(False)
        # Hide axes ticks
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])

