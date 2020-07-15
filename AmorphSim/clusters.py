from AmorphSim.real_sim import Cluster
from diffpy.structure import Atom
from AmorphSim.utils.vector_utils import get_ico_edges, get_ico_faces
from itertools import combinations
import numpy as np


class Icosahedron(Cluster):
    def __init__(self, central_atom=None, outer_atoms=None, size=1, anti=False):
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
        v = [v / np.linalg.norm(v) if np.linalg.norm(v)!= 0 else v for v in vertices]
        all_atoms = v
        if size > 1:
            if not anti:
                faces, vectors = get_ico_faces(v[1:])
                for num_shell in range(2, size+1):
                    shell = v[1:] * num_shell
                    shell2 = []
                    print("Numshells", num_shell)
                    for face in faces:
                        unique_combos = np.unique(list(combinations(list(face) * (num_shell), num_shell)), axis=0)
                        for c in unique_combos:
                            new_atom = np.sum([shell[atom - 1]*num_shell for atom in c], axis=0)/num_shell
                            shell2.append(new_atom)
                    shell=np.append(shell, shell2, axis=0)
                    all_atoms = np.append(all_atoms, shell, axis=0)
            else:
                faces, vectors = get_ico_faces(v[1:])
                for num_shell in range(2, size+1):
                    shell = np.multiply(v[1:], num_shell)
                    print(num_shell)
                    shell2 = []
                    print("Numshells", num_shell)
                    for face in faces:
                        unique_combos = np.unique(list(combinations(list(face)*(num_shell-1), num_shell+1)),axis=0)
                        for c in unique_combos:
                            new_atom = np.sum([shell[atom - 1] for atom in c], axis=0)/(num_shell+1)
                            shell2.append(new_atom)
                    shell = np.append(shell, shell2, axis=0)
                    all_atoms = np.append(all_atoms, shell, axis=0)
        all_atoms = np.unique(all_atoms,axis=0)
        print("Total Atoms : ", len(all_atoms))
        for a in all_atoms:
            self.append(Atom(xyz=a, atype=central_atom))


class FCC(Cluster):
    def __init__(self, atom1=None, atom2=None, lattice_parameter=2.54, size=1):
        """
        Parameters
        ______________
        atom_list: str or list
            The element or elements of the cluster
        lattice_parameter: float
            The lattice parameter of the cluster in angstroms
        radius: float
            The radius of the cluster in angstroms
        """
        super().__init__()
        x_pos, y_pos, z_pos = np.mgrid[-size:size + 1, -size:size + 1, -size:size + 1]
        x_pos = x_pos*lattice_parameter
        y_pos =y_pos*lattice_parameter
        z_pos = z_pos*lattice_parameter
        dis = [[0, 0, 0], [.5, .5, 0], [0, .5, .5], [.5, 0, .5]]
        for x, y, z in zip(x_pos.flatten(), y_pos.flatten(), z_pos.flatten()):
            for d in dis:
                if d == [0, 0, 0] and np.linalg.norm([x, y, z]) < size:
                    self.append(Atom(xyz=[x, y, z], atype=atom1))
                elif np.linalg.norm([x+d[0]*lattice_parameter,
                                     y+d[1]*lattice_parameter,
                                     z+d[2]*lattice_parameter]) < size:
                    self.append(Atom(xyz=[x+d[0]*lattice_parameter,
                                          y+d[1]*lattice_parameter,
                                          z+d[2]*lattice_parameter],
                                     atype=atom1))
                else:
                    pass

    def get_6_fold_axis(self, inplace=True):
        self.rotate_from_vectors(vector1=[0, 0, 1], vector2=[1, 1, 1], inplace=inplace)

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
        x_pos, y_pos, z_pos = np.mgrid[-size:(size + 1), -size:(size + 1), -size:(size + 1)]
        for x, y, z in zip(x_pos.flatten(), y_pos.flatten(), z_pos.flatten()):
            if np.linalg.norm([x, y, z]) < size:
                a1 = Atom(xyz=[x, y, z], atype=atom1)
                self.append(a1)
            if np.linalg.norm([x+.5, y+.5, z+.5]) < size:
                a2 = Atom(xyz=[x+0.5, y+0.5, z+0.5], atype=atom2)
                self.append(a2)