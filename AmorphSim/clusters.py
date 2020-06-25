from AmorphSim.real_sim import Cluster
from diffpy.structure import Atom
from AmorphSim.utils.vector_utils import get_ico_edges,get_ico_faces
import numpy as np

class Icosohedron(Cluster):
    def __init__(self, central_atom=None, outer_atoms=None):
        """
        Parameters
        ______________
        central_atom: str
            The element of the center atom
        outer_atoms: str or list
            The element of list of elements for the outer atoms of the icosohedron
        """
        super().__init__()
        golden_ratio = .5*(1+np.sqrt(5))
        vertices = [[0, 0, 0],
                    [0, 1, golden_ratio], [0, -1, golden_ratio], [0, 1, -golden_ratio], [0, -1, -golden_ratio],
                    [1, golden_ratio, 0], [-1, golden_ratio, 0], [1, -golden_ratio, 0], [-1, -golden_ratio, 0],
                    [golden_ratio, 0, 1], [-golden_ratio, 0, 1], [golden_ratio, 0, -1], [-golden_ratio, 0, -1]]
        center_atom = Atom(xyz=vertices[0], atype=central_atom)
        self.append(central_atom)
        for v in vertices[1:]:
            new_v = v / np.linalg.norm(v)
            self.append(Atom(xyz=new_v, atype=outer_atoms))


class Mackay(Icosohedron):
    def __init__(self, central_atom=None, outer_atoms=None, num_shells=2):
        """
        Parameters
        ______________
        central_atom: str
            The element of the center atom
        outer_atoms: str or list
            The element of list of elements for the outer atoms of the icosohedron
        """
        super().__init__(central_atom=central_atom,outer_atoms=outer_atoms)
        new_pos = []
        shell2 = self.xyz[1:]*2
        shell2 = np.append(shell2, get_ico_edges(shell2),axis=0)
        for atom in shell2:
            self.append(Atom(atype=outer_atoms, xyz=atom))

class AntiMackay(Icosohedron):
    def __init__(self, central_atom=None, outer_atoms=None, num_atoms=55):
        """
        Parameters
        ______________
        central_atom: str
            The element of the center atom
        outer_atoms: str or list
            The element of list of elements for the outer atoms of the icosohedron
        """
        super().__init__(central_atom=central_atom, outer_atoms=outer_atoms)
        new_pos = []
        shell2 = self.xyz[1:] * 2
        shell2 = np.append(shell2, get_ico_faces(shell2), axis=0)
        for atom in shell2:
            self.append(Atom(atype=outer_atoms, xyz=atom))

class FCC(Cluster):
    def __init__(self, atom1=None, atom2=None, lattice_parameter=1, size=1):
        """
        Parameters
        ______________
        atom_list: str or list
            The element or elements of the cluster
        lattice_parameter: float
            The lattice parameter of the cluster
        size: float
            The size of the cluster
        """
        super().__init__()
        x_pos, y_pos, z_pos = np.mgrid[-size:size + 1, -size:size + 1, -size:size + 1]
        for x, y, z in zip(x_pos.flatten(), y_pos.flatten(), z_pos.flatten()):
            a1 = Atom(xyz=[x, y, z], atype=atom1)
            a2 = Atom(xyz=[x + 0.5, y + 0.5, z], atype=atom2)
            a3 = Atom(xyz=[x + 0.5, y, z + 0.5], atype=atom2)
            a4 = Atom(xyz=[x, y + 0.5, z + 0.5], atype=atom2)
            self.append(a1)
            self.append(a2)
            self.append(a3)
            self.append(a4)

class BCC(Cluster):
    def __init__(self, atom1=None, atom2=None, lattice_parameter=1, size=1):
        """
        Parameters
        ______________
        atom_list: str or list
            The element or elements of the cluster
        lattice_parameter: float
            The lattice parameter of the cluster
        size: float
            The size of the cluster
        """
        super().__init__()
        x_pos, y_pos, z_pos = np.mgrid[-size/2:(size + 1)/2, -size/2:(size + 1)/2, -size/2:(size + 1)/2]
        for x, y, z in zip(x_pos.flatten(),y_pos.flatten(), z_pos.flatten()):
            a1 = Atom(xyz=[x, y, z], atype=atom1)
            a2 = Atom(xyz=[x+0.5, y+0.5, z+0.5], atype=atom2)
            self.append(a1)
            self.append(a2)

class Cubic(Cluster):
    def __init__(self, atom_list=None, lattice_parameter=1, size=1):
        """
        Parameters
        ______________
        atom_list: str or list
            The element or elements of the cluster
        lattice_parameter: float
            The lattice parameter of the cluster
        size: float
            The size of the cluster
        """
        pass

class Hexagonal(Cluster):
    def __init__(self, atom_list=None, lattice_parameter=1, size=1):
        """
        Parameters
        ______________
        atom_list: str or list
            The element or elements of the cluster
        lattice_parameter: float
            The lattice parameter of the cluster
        size: float
            The size of the cluster
        """
        pass