from AmorphSim.real_sim import Cluster
from diffpy.structure import Atom
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
            self.append(Atom(xyz=v, atype=outer_atoms))


class Super_Icosohedron(Cluster):
    def __init__(self, central_atom=None, outer_atoms=None):
        """
        Parameters
        ______________
        central_atom: str
            The element of the center atom
        outer_atoms: str or list
            The element of list of elements for the outer atoms of the icosohedron
        """
        pass

class FCC(Cluster):
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

class BCC(Cluster):
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