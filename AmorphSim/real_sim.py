from builtins import len

import numpy as np
from diffpy.structure import Atom, Structure
from mpl_toolkits.mplot3d import Axes3D
from AmorphSim.utils.rotation_utils import _rand_2d_rotation_matrix,_rand_3d_rotation_matrix, _get_points_on_sphere
import os
from AmorphSim.utils.vector_utils import rotation_matrix_from_vectors
import matplotlib.pyplot as plt
element_dict = {"H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8,"F": 9, "Ne": 10,
                "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18, "K": 19,
                "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28,
                "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36, "Rb": 37,
                "Sr": 38, "Y": 39, "Zr": 40}


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
        return "<Cube of " + len(self.clusters) + " clusters"

    def add_cluster(self, cluster, location):
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
        ax.scatter(coordinates[:, 0], coordinates[:,1], coordinates[:,2],
                   marker='o', s=300)

    def random_rot_xyz(self, num=100, folder="RandCluster", scale=2, offset=[100,100,20]):
        if not os.path.exists(folder):
            os.makedirs(folder)
        for i in range(num):
            M = _rand_3d_rotation_matrix()
            new = self.rotate_from_matrix(matrix=M, inplace=False)
            new.get_xyz(file=folder+"/"+str(i)+".xyz",scale=scale, offset=offset)

    def get_xyz(self, file=None, scale=5, offset=[100, 100, 20]):
        positions = self.xyz_cartn
        atom_list = self.tolist()
        newstr = "Simulating XYZ \n     200.0 200.0 40.0\n"
        for a in atom_list:
            newstr = newstr+(str(element_dict[a.element])+ " " +
                             str((a.xyz[0]*scale)+offset[0]) + " " +
                             str((a.xyz[1]*scale)+offset[1])+ " " +
                             str((a.xyz[2]*scale)+offset[2])+ " " + str(1)+ " " + str(0.5)+"\n")
        if file is None:
            return newstr
        else:
            with open(file, "+w") as f:
                f.write(newstr)

    def rotate_from_matrix(self, matrix, inplace=False):
        if inplace:
            for i in range(len(self)):
                self[i].xyz = np.dot(matrix, self[i].xyz)
            return
        else:
            new = self
            for i in range(len(new)):
                new[i].xyz = np.dot(matrix, new[i].xyz)
            return new

    def all_direction_xyz(self,npt=1000, folder="alldirections", scale=2,offset=[100,100,20]):
        vectors = _get_points_on_sphere(npt=npt)
        if not os.path.exists(folder):
            os.makedirs(folder)
        for v in vectors:
            rotated =self.rotate_from_vectors(vector1=[1,0,0],
                                              vector2=v,
                                              inplace=False)
            print(rotated)
            rotated.get_xyz(file=folder + "/" + str(v)+ ".xyz",
                            scale=scale,
                            offset=offset)


    def rotate_from_vectors(self, vector1, vector2, inplace=False):
        mat = rotation_matrix_from_vectors(vec1=vector1, vec2=vector2)
        if inplace:
            for i in range(len(self)):
                self[i].xyz = np.dot(mat, self[i].xyz)
            return None
        else:
            new = self
            for i in range(len(new)):
                new[i].xyz = np.dot(mat, new[i].xyz)
            return new

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
