import numpy as np
from diffpy.structure import Atom
from AmorphSim.utils.vector_utils import get_ico_edges, get_ico_faces
from itertools import combinations


def _create_ico(central_atom="Cu",
                outer_atoms="Cu",
                shell_distance=2.54,
                num_shells=10,
                anti=False):
    """Sets the initial positions for an icosahedron-like cluster.
      Uses the golden ratio to create a set of evenly
      spaced atom around some central point.

    Parameters
    _______________
    central_atom: str
        The central atom of the icosahedron or Mackay cluster
    outer_atoms: str
        The outer atoms of the icosahedron of Mackay cluster
    shell_distance: float
        The distance between shells of the Mackay cluster
    num_shells: int
        The number of shells in the cluster
    anti: bool
        If the cluster is Anti-Mackay cluster

    Returns
    ------------
    atom_list: list
        A list of atoms in the cluster
    """
    atom_list = []
    golden_ratio = .5 * (1 + np.sqrt(5))
    vertices = [[0, 0, 0],
                [0, 1, golden_ratio], [0, -1, golden_ratio], [0, 1, -golden_ratio], [0, -1, -golden_ratio],
                [1, golden_ratio, 0], [-1, golden_ratio, 0], [1, -golden_ratio, 0], [-1, -golden_ratio, 0],
                [golden_ratio, 0, 1], [-golden_ratio, 0, 1], [golden_ratio, 0, -1], [-golden_ratio, 0, -1]]
    v = [(v / np.linalg.norm(v)) * shell_distance if np.linalg.norm(v) != 0 else v for v in vertices]
    all_atoms = v
    if num_shells > 1:
        if not anti:
            faces, vectors = get_ico_faces(v[1:], shell_dist=shell_distance)
            for n_shell in range(2, num_shells + 1):
                shell = v[1:] * n_shell
                shell2 = []
                print("Numshells", n_shell)
                for face in faces:
                    unique_combos = np.unique(list(combinations(list(face) * (n_shell), n_shell)), axis=0)
                    for c in unique_combos:
                        new_atom = np.sum([shell[atom - 1] * n_shell for atom in c], axis=0) / n_shell
                        shell2.append(new_atom)
                shell = np.append(shell, shell2, axis=0)
                all_atoms = np.append(all_atoms, shell, axis=0)
        else:
            faces, vectors = get_ico_faces(v[1:], shell_dist=shell_distance)
            for n_shell in range(2, num_shells + 1):
                shell = np.multiply(v[1:], n_shell)
                print(n_shell)
                shell2 = []
                print("Numshells", n_shell)
                for face in faces:
                    unique_combos = np.unique(list(combinations(list(face) * (n_shell - 1), n_shell + 1)), axis=0)
                    for c in unique_combos:
                        new_atom = np.sum([shell[atom - 1] for atom in c], axis=0) / (n_shell + 1)
                        shell2.append(new_atom)
                shell = np.append(shell, shell2, axis=0)
                all_atoms = np.append(all_atoms, shell, axis=0)
    all_atoms = np.unique(all_atoms, axis=0)
    print("Total Atoms : ", len(all_atoms))
    for a in all_atoms:
        atom_list.append(Atom(xyz=a, atype=central_atom))
    return atom_list


def _create_fcc(atom1=None, atom2=None, lattice_parameter=2.54, size=10):
    """Sets the initial positions for an FCC-like cluster.
      Uses the golden ratio to create a set of evenly
      spaced atom around some central point.

    Parameters
    _______________
    atom1: str
        The corner atom of the cluster
    atom2: str
        The edge atoms of the cluster
    lattice_parameter: float
        The lattice parameter of the unit cell
    size: float
        The size of the cluster in angstroms

    Returns
    ------------
    atom_list: list
        A list of atoms in the cluster
    """
    atom_list = []
    x_pos, y_pos, z_pos = np.mgrid[-size:size + 1, -size:size + 1, -size:size + 1]
    x_pos = x_pos * lattice_parameter
    y_pos = y_pos * lattice_parameter
    z_pos = z_pos * lattice_parameter
    dis = [[0, 0, 0], [.5, .5, 0], [0, .5, .5], [.5, 0, .5]]
    for x, y, z in zip(x_pos.flatten(), y_pos.flatten(), z_pos.flatten()):
        for d in dis:
            if d == [0, 0, 0] and np.linalg.norm([x, y, z]) < size:
                atom_list.append(Atom(xyz=[x, y, z], atype=atom1))
            elif np.linalg.norm([x + d[0] * lattice_parameter,
                                 y + d[1] * lattice_parameter,
                                 z + d[2] * lattice_parameter]) < size:
                atom_list.append(Atom(xyz=[x + d[0] * lattice_parameter,
                                      y + d[1] * lattice_parameter,
                                      z + d[2] * lattice_parameter],
                                 atype=atom2))
            else:
                pass
    return atom_list


def _create_bcc(atom1=None, atom2=None, lattice_parameter=2.54, size=10):
    """Sets the initial positions for an FCC-like cluster.
      Uses the golden ratio to create a set of evenly
      spaced atom around some central point.

    Parameters
    _______________
    atom1: str
        The corner atom of the cluster
    atom2: str
        The edge atoms of the cluster
    lattice_parameter: float
        The lattice parameter of the unit cell
    size: float
        The size of the cluster in angstroms

    Returns
    ------------
    atom_list: list
        A list of atoms in the cluster
    """
    atom_list = []
    int_size = np.ceil(size/lattice_parameter)
    print(int_size)
    x_pos, y_pos, z_pos = np.mgrid[-int_size:(int_size + 1),
                                   -int_size:(int_size + 1),
                                   -int_size:(int_size + 1)]
    x_pos = x_pos * lattice_parameter
    y_pos = y_pos * lattice_parameter
    z_pos = z_pos * lattice_parameter
    for x, y, z in zip(x_pos.flatten(), y_pos.flatten(), z_pos.flatten()):
        if np.linalg.norm([x, y, z]) < size:
            a1 = Atom(xyz=[x, y, z], atype=atom1)
            atom_list.append(a1)
        if np.linalg.norm([x + .5*lattice_parameter,
                           y + .5*lattice_parameter,
                           z + .5*lattice_parameter]) < size:
            a2 = Atom(xyz=[x + 0.5*lattice_parameter,
                           y + 0.5*lattice_parameter,
                           z + 0.5*lattice_parameter],
                      atype=atom2)
            atom_list.append(a2)
    return atom_list
