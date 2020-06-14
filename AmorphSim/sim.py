import numpy as np
#import hyperspy.api as hs
#from hyperspy._signals.signal2d import Signal2D

from AmorphSim.utils.rotation_utils import _get_rotation_matrix, _get_random_2d_rot, _get_random_3d_rot
from AmorphSim.utils.simulation_utils import _get_speckle_size, _get_wavelength, _shape_function, _get_speckle_intensity
from AmorphSim.utils.vector_utils import rotation_matrix_from_vectors, build_ico
from skimage.draw import circle
from skimage.filters import gaussian
from numpy.random import random, choice
import matplotlib.pyplot as plt
from matplotlib.patches import Circle


class SimulationCube(object):
    """Defines a simulation cube of dimensions x, y, z in nm.  This allows you to create some simulation of the cube
    based on kinematic diffraction"""
    def __init__(self, dimensions=(20, 20, 20)):
        """Initializes the simulation cube
        Parameters
        --------------
        dimensions: tuple
            The dimensions of the cube to simulate.
        """
        self.clusters = []
        self.dimensions = np.array(dimensions)

    def __str__(self):
        return ("<Cube of " + str(self.num_clusters) + " clusters [" +
               str(self.dimensions[0])+" x "+str(self.dimensions[1])
               +" x "+str(self.dimensions[2])+"nm]>")

    @property
    def num_clusters(self):
        return len(self.clusters)

    def add_random_clusters(self, num_clusters, radius_range=(.5, 1.0), k_range=(3.5, 4.5), random_rotation=True,
                            symmetry=[2, 4, 6, 10]):
        """Randomly initializes the glass with a set of random clusters.
        Parameters
        ------------
        num_clusters: int
            The number of cluster to add
        radius_range: tuple
            The range of radii to randomly choose from
        k_range: tuple
            THe range of k to randomly choose from
        random_rotation: bool
            Random rotate the cluster or not
        symmetry: list
            The list of symmetries to choose from.  Acceptable symmetries are 2,4,6 and 10
        """
        rand_r = random(num_clusters) * (radius_range[1] - radius_range[0]) + radius_range[0]
        rand_k = random(num_clusters) * (k_range[1] - k_range[0]) + k_range[0]
        rand_sym = choice(symmetry,num_clusters)
        rand_pos = np.multiply(random((num_clusters,3)), self.dimensions)
        if random_rotation:
            rot_3d = _get_random_3d_rot(num_clusters)
            rot_2d = _get_random_2d_rot(num_clusters)
        else:
            rot_2d = [np.eye(3), ]*num_clusters
            rot_3d = [np.eye(3), ]*num_clusters
        for s, r, k, p, two, three in zip(rand_sym,rand_r, rand_k,rand_pos, rot_2d, rot_3d):
            self.clusters.append(Cluster(s,r,k,p,two,three))
        return

    def add_icoso(self, num_clusters, radius_range=(.5, 1.0), k_range=(3.5, 4.5), random_rotation=True):
        """Randomly initializes the glass with a set of random clusters.
        Parameters
        ------------
        num_clusters: int
            The number of cluster to add
        radius_range: tuple
            The range of radii to randomly choose from
        k_range: tuple
            THe range of k to randomly choose from
        random_rotation: bool
            Random rotate the cluster or not
        symmetry: list
            The list of symmetries to choose from.  Acceptable symmetries are 2,4,6 and 10
        """
        rand_r = random(num_clusters) * (radius_range[1] - radius_range[0]) + radius_range[0]
        rand_k = random(num_clusters) * (k_range[1] - k_range[0]) + k_range[0]
        rand_pos = np.multiply(random((num_clusters,3)), self.dimensions)
        if random_rotation:
            rot_3d = _get_random_3d_rot(num_clusters)
            rot_2d = _get_random_2d_rot(num_clusters)
        else:
            rot_2d = [np.eye(3), ]*num_clusters
            rot_3d = [np.eye(3), ]*num_clusters
        five_vertexes, three_face, two_edge = build_ico()
        for r, k, p, two, three in zip(rand_r, rand_k,rand_pos, rot_2d, rot_3d):
            for v in five_vertexes:
                self.clusters.append(Cluster(10,r,k,p,two,three,plane_direction=v))
            for v in three_face:
                self.clusters.append(Cluster(6,r,k,p,two,three,plane_direction=v))
            for v in two_edge:
                self.clusters.append(Cluster(2, r, k, p, two, three, plane_direction=v))
        return

    def show_projection(self, size=(512,512), acceptance = None):
        """Plots the 2-d projection.  Creates a 2-D projection of the clusters in the amorphous matrix.

        Parameters:
        ------------
        acceptance: float
            The angle of acceptance.  Only clusters within this projection will be allowed.
        size: int
            The size of the image made

        """
        projection = np.zeros(size)
        scale = [self.dimensions[0]/size[0],self.dimensions[1]/size[1]]
        for cluster in self.clusters:
            if acceptance is not None:
                if abs(cluster.get_angle_between()) > acceptance:
                    next()
            inten = np.sum(cluster.get_intensity())
            r, c = circle(cluster.position[0]/scale[0],
                          cluster.position[1]/scale[1],
                          radius=cluster.radius/scale[0],
                          shape=size)
            projection[r, c] = inten + projection[r, c]
        return projection

    def plot_symmetries(self, symmetries=[2, 4, 6, 8, 10], norm=1,acceptance=None):
        """Plots the 2-d projection of the symmetries as circles

        Parameters:
        ------------
        acceptance: float
            The angle of acceptance.  Only clusters within this projection will be allowed.
        size: int
            The size of the image made

        """
        fig, ax = plt.subplots()
        ax.set_xlim(0,self.dimensions[0])
        ax.set_ylim(0,self.dimensions[1])
        colors = ["black","blue","red","green","yellow","red","orange", "purple"]
        for cluster in self.clusters:
            if acceptance is not None:
                if abs(cluster.get_angle_between()) > acceptance:
                    continue
            inten = np.mean(cluster.get_intensity())
            c = Circle((cluster.position[0],
                        cluster.position[1]),
                        radius=cluster.radius,
                        alpha=inten/norm,
                        color=colors[symmetries.index(cluster.symmetry)])
            ax.add_patch(c)
        from matplotlib.lines import Line2D
        leg = [Line2D([0], [0], marker='o', color=colors[i], label=str(sym)+" fold symmetry",
               markerfacecolor=colors[i], markersize=15) for i,sym in enumerate(symmetries)]

        ax.legend(handles=leg)
        plt.show()
        return

    def get_4d_stem(self, convergence_angle=.74, accelerating_voltage=200,
                    k_rad = 5.0, simulation_size=(50, 50, 128, 128),
                    noise = False, disorder = None, num_electrons=1000, convolve=False, beam_size=None):
        """Returns an amorphous2d object which shows the 4d STEM projection for some set of clusters along some
        illumination

        Parameters
        ------------
        convergence_angle: float
            The convergance angle for the experiment
        accelerating_voltage: float
            The accelerating voltage for the experiment in kV
        simulation_size: tuple
            The size of the image for both the reciporical space image and the real space image.

        Returns
        ------------
        dataset: Amorphus2D
            Returns a 4 dimensional dataset which represents the cube
        """
        # dataset = Signal2D(np.ones(simulation_size))
        dataset = np.ones(simulation_size)
        real_scale = simulation_size[0] / self.dimensions[0]

        for cluster in self.clusters:
            speckles, observed_intensity = cluster.get_speckles(img_size=k_rad*2,
                                                                num_pixels=simulation_size[2],
                                                                accelerating_voltage=accelerating_voltage,
                                                                conv_angle=convergence_angle,
                                                                disorder=disorder)

            rr, rc = circle(r=cluster.position[0] * real_scale,
                            c=cluster.position[1] * real_scale,
                            radius=cluster.radius * real_scale,
                            shape=(simulation_size[0], simulation_size[1]))
            for (sr,sc), inten in zip(speckles, observed_intensity):
                inner_r, outer_r = np.meshgrid(sr, rr)
                inner_c, outer_c = np.meshgrid(sc, rc)
                if noise:
                    #print(np.size(inner_c))
                    inten = np.random.poisson(inten*num_electrons,size=np.size(inner_c))
                else:
                    inten= inten* num_electrons
                dataset[(outer_c.flatten(),
                         outer_r.flatten(),
                         inner_c.flatten(),
                         inner_r.flatten())] = inten + dataset[(outer_c.flatten(),
                                                                outer_r.flatten(),
                                                                inner_c.flatten(),
                                                                inner_r.flatten())]
        dataset = Signal2D(dataset)
        if convolve:
            if beam_size is None:
                beam_size = .9/convergence_angle  # rough calculation
            from scipy.signal import convolve2d
            num_pixels = int(real_scale*beam_size)
            xx,yy= np.ogrid[-num_pixels-1:num_pixels+2,-num_pixels-1:num_pixels+2]
            kernel =1 - ((xx ** 2 + yy ** 2) ** .5 - num_pixels).clip(0,1)
            dataset = dataset.T
            dataset.map(convolve2d, in2=kernel, mode="same", inplace=True)
            dataset = dataset.T


        dataset.axes_manager.navigation_axes[0].scale = self.dimensions[0]/simulation_size[0]
        dataset.axes_manager.navigation_axes[1].scale = self.dimensions[1] / simulation_size[1]
        dataset.axes_manager.navigation_axes[0].units = "nm"
        dataset.axes_manager.navigation_axes[1].units = "nm"
        dataset.axes_manager.signal_axes[0].scale = k_rad*2 /simulation_size[2]
        dataset.axes_manager.signal_axes[1].scale = k_rad*2 / simulation_size[3]
        dataset.axes_manager.signal_axes[0].units = "$nm^-1$"
        dataset.axes_manager.signal_axes[1].units = "$nm^-1$"
        dataset.axes_manager.signal_axes[0].offset = -k_rad
        dataset.axes_manager.signal_axes[1].offset = -k_rad
        return dataset


class Cluster(object):
    def __init__(self,
                 symmetry=10,
                 radius=1,
                 k=4.0,
                 position=random(2),
                 rotation_2d=np.eye(3),
                 rotation_3d=np.eye(3),
                 plane_direction=[0,0,1]):
        """Defines a cluster with a symmetry of symmetry, a radius of radius in nm and position of position.

        Parameters:
        ----------------
        symmetry: int
            The symmetry of the cluster being simulated
        radius: float
            The radius of the cluster in nm
        position: tuple
            The position of the cluster in the simulation cube
        rotation_vector: tuple
            The vector which the cluster is rotated around
        rotation_angle: float
            The angle the cluster is rotated about
        """
        self.symmetry = symmetry
        self.radius = radius
        self.position = position
        self.rotation_2d = rotation_2d
        self.rotation_3d = rotation_3d
        self.k = k
        self.plane_direction = plane_direction
        self.beam_direction = [0,0,1]

    def get_diffraction(self, img_size=8.0,
                        num_pixels=512,
                        accelerating_voltage=200,
                        conv_angle=0.6,
                        disorder=None):
        """Takes some image size in inverse nm and then plots the resulting
        """
        sphere_radius = 1/_get_wavelength(accelerating_voltage)
        scale = (num_pixels-1)/img_size
        k_rotated=self.get_k_vectors()
        observed_intensity = [_get_speckle_intensity(k_vector=k,
                                                     ewald_sphere_rad=sphere_radius,
                                                     disorder=disorder,
                                                     cluster_rad=self.radius,
                                                     beam_direction=self.beam_direction)
                              for k in k_rotated]
        radius = _get_speckle_size(accelerating_voltage, conv_angle)*scale
        circles = [circle(int(k1[0] * scale + num_pixels/2), int(k1[1] * scale + num_pixels/2),
                          radius=radius) for k1 in k_rotated]
        image = np.ones(shape=(num_pixels,num_pixels))
        for (r, c), i in zip(circles, observed_intensity):
            image[r, c] = i + image[r, c]
        image = gaussian(image=image, sigma=2)
        return image

    def get_k_vectors(self):
        angle = (2 * np.pi) / self.symmetry  # angle between speckles on the pattern
        k = [[np.cos(angle * i) * self.k, np.sin(angle * i) * self.k, 0] for i in
             range(self.symmetry)]  # vectors for the speckles perp to BA
        #print(self.plane_direction)
        if not np.allclose(self.plane_direction, [0, 0, 1]) and not np.allclose(self.plane_direction, [0, 0, -1]) :
            rot = rotation_matrix_from_vectors([0,0,1],self.plane_direction)
            k = [np.dot(rot,v) for v in k]
        k_rotated2d = [np.dot(self.rotation_2d, speckle) for speckle in k]
        k_rotated = [np.dot(self.rotation_3d, speckle) for speckle in k_rotated2d]
        return k_rotated

    def get_speckles(self, img_size=10.0,
                     num_pixels=128,
                     accelerating_voltage=200,
                     conv_angle=0.6,
                     disorder=None):
        """
        This function returns the diffraction speckles as circles as defined by the
        skimage.draw.Circle class. Each speckle also has some intensity associated with it.

        Parameters
        -----------
        img_size:
            The real size of the image that is being projected onto. (In inverse nm)
        num_pixels:
            The pixelated size of the image being projected onto.
        accelerating_voltage:
            The accelerating voltage for getting the Ewald sphere radius
        conv_angle:
            The convergence angle in mrad for beam.
        """
        sphere_radius = 1/_get_wavelength(accelerating_voltage)
        scale = (num_pixels-1)/img_size
        k_rotated = self.get_k_vectors()
        observed_intensity = [_get_speckle_intensity(k_vector=k,
                                                     ewald_sphere_rad=sphere_radius,
                                                     disorder=disorder,
                                                     cluster_rad=self.radius,
                                                     beam_direction=self.beam_direction)
                              for k in k_rotated]
        radius = _get_speckle_size(accelerating_voltage, conv_angle)*scale
        speckles = [circle(int(k1[0] * scale + num_pixels/2), int(k1[1] * scale + num_pixels/2),
                           radius=radius, shape=(num_pixels,num_pixels)) for k1 in k_rotated]
        return speckles, observed_intensity

    def get_intensity(self,accelerating_voltage=200, disorder=None):
        """Takes some image size in inverse nm and then plots the resulting
        """
        sphere_radius = 1 / _get_wavelength(accelerating_voltage)
        k_rotated = self.get_k_vectors()
        observed_intensity = [_get_speckle_intensity(k_vector=k,
                                                     ewald_sphere_rad=sphere_radius,
                                                     disorder=disorder,
                                                     cluster_rad=self.radius,
                                                     beam_direction=self.beam_direction)
                              for k in k_rotated]
        return observed_intensity

    def get_angle_between(self):
        return np.arccos((np.trace(self.rotation_3d)-1)/2)
