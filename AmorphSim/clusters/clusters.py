from AmorphSim.real_sim import Cluster
from AmorphSim.utils.structure_generators import _create_ico, _create_fcc,_create_bcc
from AmorphSim.draw.draw_3d import Icosohedron_3d,FCC_3d


class Icosahedron(Cluster):
    def __init__(self,
                 central_atom="Cu",
                 outer_atoms="Cu",
                 shell_distance=2.54,
                 num_shells=1,
                 anti=False,
                 position=[0, 0, 0]):
        """
        Parameters
        ______________
        central_atom: str
            The element of the center atom
        outer_atoms: str or list
            The element of list of elements for the outer atoms of the icosahedron
        shell_distance: float
            The distance from one atom to each successive shell in the cluster.
        num_shells: int
            The number of shells in the icosahedron
        anti: boolean
            If the cluster should be an anti mackay cluster instead of a mackay cluster
        """
        super().__init__(position=position)
        self.initial_atoms = _create_ico(central_atom=central_atom,
                                         outer_atoms=outer_atoms,
                                         shell_distance=shell_distance,
                                         num_shells=num_shells,
                                         anti=anti)
        for i in self.initial_atoms:
            self.append(i)
        self.radius = num_shells*shell_distance

    def get_5_fold_axis(self, beam_direction=[0, 0, 1], reinitialize=True, inplace=True):
        """rotate the 5fold axes so that it is perpendicular to the beam direction [0,0,1]

        Parameters
        ------------
        beam_direction: list
            The beam direction to align the axis to
        reinitialize: bool
            Resets the atom positions so that the axes is aligned properly
        inplace: bool
            Should a copy be made or the cluster be directly operated on
        """
        if reinitialize:
            self.reset_atoms()
        return self.rotate_from_vectors(vector1=beam_direction,
                                        vector2=[0.0, -0.52573111, 0.85065081],
                                        inplace=inplace)

    def get_3_fold_axis(self, beam_direction=[0, 0, 1], reinitialize=True, inplace=True):
        """rotate to 3-fold axes so that it is perpendicular to the beam direction [0,0,1]

        Parameters
        ------------
        beam_direction: list
            The beam direction to align the axis to
        reinitialize: bool
            Resets the atom positions so that the axes is aligned properly
        inplace: bool
            Should a copy be made or the cluster be directly operated on
        """
        if reinitialize:
            self.reset_atoms()
        return self.rotate_from_vectors(vector1=beam_direction,
                                        vector2=[-0.45879397, -0.45879397, 0.45879397],
                                        inplace=inplace)

    def get_2_fold_axis(self, beam_direction=[0, 0, 1], reinitialize=True, inplace=True):
        """rotate to 3-fold axes so that it is perpendicular to the beam direction [0,0,1]

        Parameters
        ------------
        beam_direction: list
            The beam direction to align the axis to
        reinitialize: bool
            Resets the atom positions so that the axes is aligned properly
        inplace: bool
            Should a copy be made or the cluster be directly operated on
        """
        if reinitialize:
            self.reset_atoms()
        if inplace:
            return None
        else:
            return self

    def draw(self):
        print("Radius:", self.radius)
        return Icosohedron_3d(radius=self.radius, location=self.position)


class FCC(Cluster):
    def __init__(self,
                 atom1="Cu",
                 atom2="Cu",
                 lattice_parameter=2.54,
                 radius=3,
                 position=[0, 0, 0]):
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
        super().__init__(position=position)
        self.initial_atoms = _create_fcc(atom1=atom1,
                                         atom2=atom2,
                                         lattice_parameter=lattice_parameter,
                                         size=radius)
        for i in self.initial_atoms:
            self.append(i)
        self.radius = radius

    def get_6_fold_axis(self, beam_direction=[0, 0, 1], reinitialize=True, inplace=True):
        """rotate to 6-fold axes so that it is perpendicular to the beam direction [0,0,1]

        Parameters
        ------------
        beam_direction: list
            The beam direction to align the axis to
        reinitialize: bool
            Resets the atom positions so that the axes is aligned properly
        inplace: bool
            Should a copy be made or the cluster be directly operated on
        """
        if reinitialize:
            self.reset_atoms()
        return self.rotate_from_vectors(vector1=beam_direction,
                                        vector2=[1, 1, 1],
                                        inplace=inplace)

    def get_4_fold_axis(self, beam_direction=[0, 0, 1], reinitialize=True, inplace=True):
        """rotate to 4-fold axes so that it is perpendicular to the beam direction [0,0,1]

        Parameters
        ------------
        beam_direction: list
            The beam direction to align the axis to
        reinitialize: bool
            Resets the atom positions so that the axes is aligned properly
        inplace: bool
            Should a copy be made or the cluster be directly operated on
        """
        if reinitialize:
            self.reset_atoms()
        if inplace:
            return None
        else:
            return self

    def draw(self):
        return FCC_3d(radius=self.radius, location=self.position)


class BCC(Cluster):
    def __init__(self, atom1="Fe", atom2="Fe", lattice_parameter=1, radius=1):
        """
        Parameters
        ______________
        atom_list: str or list
            The element or elements of the cluster
        lattice_parameter: float
            The lattice parameter of the cluster
        radius: float
            The size of the cluster
        """
        super().__init__()
        self.initial_atoms = _create_bcc(atom1=atom1,
                                         atom2=atom2,
                                         lattice_parameter=lattice_parameter,
                                         size=radius)
        for a in self.initial_atoms:
            self.append(a)

    def get_6_fold_axis(self, beam_direction=[0, 0, 1], reinitialize=True, inplace=True):
        """rotate to 6-fold axes so that it is perpendicular to the beam direction [0,0,1]

        Parameters
        ------------
        beam_direction: list
            The beam direction to align the axis to
        reinitialize: bool
            Resets the atom positions so that the axes is aligned properly
        inplace: bool
            Should a copy be made or the cluster be directly operated on
        """
        if reinitialize:
            self.reset_atoms()
        return self.rotate_from_vectors(vector1=beam_direction,
                                        vector2=[1, 1, 1],
                                        inplace=inplace)

    def get_4_fold_axis(self, beam_direction=[0, 0, 1], reinitialize=True, inplace=True):
        """rotate to 4-fold axes so that it is perpendicular to the beam direction [0,0,1]

        Parameters
        ------------
        beam_direction: list
            The beam direction to align the axis to
        reinitialize: bool
            Resets the atom positions so that the axes is aligned properly
        inplace: bool
            Should a copy be made or the cluster be directly operated on
        """
        if reinitialize:
            self.reset_atoms()
        if inplace:
            return None
        else:
            return self


class BCC(Cluster):
    def __init__(self, atom1="Fe", atom2="Fe", lattice_parameter=2.856, radius=3):
        """
        Parameters
        ______________
        atom_list: str or list
            The element or elements of the cluster
        lattice_parameter: float
            The lattice parameter of the cluster
        radius: float
            The size of the cluster
        """
        super().__init__()
        self.initial_atoms = _create_bcc(atom1=atom1,
                                         atom2=atom2,
                                         lattice_parameter=lattice_parameter,
                                         size=radius)
        for a in self.initial_atoms:
            self.append(a)

    def get_6_fold_axis(self, beam_direction=[0, 0, 1], reinitialize=True, inplace=True):
        """rotate to 6-fold axes so that it is perpendicular to the beam direction [0,0,1]

        Parameters
        ------------
        beam_direction: list
            The beam direction to align the axis to
        reinitialize: bool
            Resets the atom positions so that the axes is aligned properly
        inplace: bool
            Should a copy be made or the cluster be directly operated on
        """
        if reinitialize:
            self.reset_atoms()
        return self.rotate_from_vectors(vector1=beam_direction,
                                        vector2=[1, 1, 1],
                                        inplace=inplace)

    def get_4_fold_axis(self, beam_direction=[0, 0, 1], reinitialize=True, inplace=True):
        """rotate to 4-fold axes so that it is perpendicular to the beam direction [0,0,1]

        Parameters
        ------------
        beam_direction: list
            The beam direction to align the axis to
        reinitialize: bool
            Resets the atom positions so that the axes is aligned properly
        inplace: bool
            Should a copy be made or the cluster be directly operated on
        """
        if reinitialize:
            self.reset_atoms()
        if inplace:
            return None
        else:
            return self