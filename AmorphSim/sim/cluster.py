from diffpy.structure import Atom, Structure,Lattice
from diffsims.generators import DiffractionGenerator
import matplotlib.pyplot as plt


element_dict = {"H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8,"F": 9, "Ne": 10,
                "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18, "K": 19,
                "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28,
                "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36, "Rb": 37,
                "Sr": 38, "Y": 39, "Zr": 40}

class Cluster(Structure):
    """Each Cluster extends the Structure class
     giving it access to the Strcuture class and
     the diffraction simulation capabilities of diffsims
    """
    def __init__(self, position=[0,0,0]):
        self.initial_atoms = []
        self.position = position

    def reset_atoms(self):
        for i, a in enumerate(self.initial_atoms):
            self[i] = a

    def plot(self, save=False, rotate=True):
        """Plots the atoms of some structure in 3-D.  The atoms are shown as spheres based on their atomic
        size.
        """
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        coordinates = self.xyz_cartn
        print(self.element)

        if not save:
            ax.scatter(coordinates[:, 0], coordinates[:, 1], coordinates[:, 2],
                       marker='o', s=300)
            ax.set_axis_off()
            if rotate:
                for angle in range(0, 360):
                    ax.view_init(30, angle)
                    plt.draw()
                    plt.pause(.001)

        if save:
            import matplotlib.animation as animation
            from matplotlib import rcParams
            def init_function():
                ax.scatter(coordinates[:, 0], coordinates[:, 1], coordinates[:, 2],
                           marker='o', s=300)
                ax.set_axis_off()
                return fig,
            def animate(i):
                # azimuth angle : 0 deg to 360 deg
                ax.view_init(elev=10, azim=i * 4)
                return fig,

            # create animation using the animate() function with no repeat
            # Animate
            ani = animation.FuncAnimation(fig, animate, init_func=init_function,
                                          frames=90, interval=50, blit=True)
            fn = 'rotate_azimuth_angle_3d_surf'
            ani.save(fn + '.mp4', writer='ffmpeg', fps=1000 / 50)
            ani.save(fn + '.gif', writer='imagemagick', fps=1000 / 50)

    def random_rot_xyz(self, num=100, folder="RandCluster", scale=2, offset=[100, 100, 20]):
        if not os.path.exists(folder):
            os.makedirs(folder)
        for i in range(num):
            M = _rand_3d_rotation_matrix()
            new = self.rotate_from_matrix(matrix=M, inplace=False)
            new.get_xyz(file=folder+"/"+str(i)+".xyz",scale=scale, offset=offset)

    def get_xyz(self, file=None, offset=[100, 100, 20], disorder=[0.05, 0.1, 0.15, 0.20, 0.25, 0.3], fp=5):
        """Get a series of xyz positions as a function of list of disorder parameters
        file: str
            The name of the file to be created without .xyz at the end of the filename
        offset: vector
            The offset for the simulation.  This depends on the simulation
            size and shifts the atom positions by the offset
        disorder: list

        """

        string_list =[]
        for d in disorder:
            newstr = "Simulating XYZ \n     200.0 200.0 40.0\n"
            for nf in range(fp):
                new = self.add_disorder(sigma=d, inplace=False)
                for a in new:
                    newstr = newstr+(str(element_dict[a.element]) + " " +
                                     str(a.x+offset[0]) + " " +
                                     str(a.y+offset[1]) + " " +
                                     str(a.z+offset[2]) + " " +
                                     str(a.occupancy) + " " +
                                     str(0.0)+"\n")
            string_list.append(newstr)
        if file is None:
            return string_list
        else:
            for d, newstr in zip(disorder, string_list):
                for nf in range(fp):
                    with open(file+"d:"+str(d)+"fp:"+str(nf)+".xyz", "+w") as f:
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

    def all_direction_xyz(self, npt=1000, folder="alldirections", scale=2, offset=[100, 100, 20], disorder=1.0):
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
                            offset=offset,
                            disorder=disorder)
        return

    def add_disorder(self, sigma=0.1, inplace=True):
        if inplace:
            for atom in self:
                dx, dy, dz = np.random.randn(3)*sigma
                atom.x = atom.x + dx
                atom.y = atom.y + dy
                atom.z = atom.z + dz
        else:
            new = self.copy()
            for atom in new:
                dx, dy, dz = np.random.randn(3)*sigma
                atom.x = atom.x + dx
                atom.y = atom.y + dy
                atom.z = atom.z + dz
            return new

    def __copy__(self, target=None):
        '''Create a deep copy of this instance.
        target   -- optional target instance for copying, useful for
                    copying a derived class.  Defaults to new instance
                    of the same type as self.
        Return a duplicate instance of this object.
        '''
        if target is None:
            target = Cluster()
        elif target is self:
            return target
        # copy attributes as appropriate:
        target.title = self.title
        target.lattice = Lattice(self.lattice)
        target.pdffit = copymod.deepcopy(self.pdffit)
        # copy all atoms to the target
        target[:] = self
        return target

    def rotate_from_vectors(self, vector1, vector2, inplace=False):
        mat = rotation_matrix_from_vectors(vec1=vector1, vec2=vector2)
        if inplace:
            for i in range(len(self)):
                self[i].xyz = np.dot(mat, self[i].xyz)
            return None
        else:
            new = self.__copy__()
            for i in range(len(new)):
                new[i].xyz = np.dot(mat, new[i].xyz)
            return new

    def plot_rotate(self, start=[1,0,0], end=[0,1,0],npts=90):
        angle = _get_angle_between(start,end)
        angles = np.linspace(0, angle,npts)
        vectors = [[np.cos(a), np.sin(a), 0] for a in angles]

    def plot_reciprocal(self, resolution=100):
        """Plots the reciprocal space  atoms of some structure in 3-D. Gives a voxel representation of
        the space.
        """

        pass

    def get_reciprocal_space(self, accelerating_voltage, resolution=100):
        """Returns the three dimensional reciprocal space from some cluster where some cut along the
        reciprocal space represents the kinetic diffraction along some direction. Creates a resolution^3
        space
        """
        gen = DiffractionGenerator(accelerating_voltage)

        pass
    def draw(self):
        return