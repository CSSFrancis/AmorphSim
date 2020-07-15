import os
from AmorphSim.real_sim import Cluster
import numpy as np
from AmorphSim.clusters import FCC, BCC, Icosahedron
from AmorphSim.utils.rotation_utils import _get_points_on_sphere


def simulate_all(atom="Cu",
                 folder="all",
                 cluster="Icosahedron",
                 npt_rot=1000,
                 sizes=[1, 2, 3],
                 disorder=[0.2, 0.4, 0.6, 0.8, 1.0]):
    if not os.path.exists(folder):
        os.makedirs(folder)
    vectors = _get_points_on_sphere(npt=npt_rot)
    if cluster == "Icosahedron":
        for s in sizes:
            i = Icosahedron(central_atom="Cu", size=s)
            for v in vectors:
                for d in disorder:
                    i.rotate_from_vectors([1, 0, 0], v, inplace=False).get_xyz(file=folder+
                                                                   "/size:"+str(s)+"disorder:"+
                                                                   str(d)+"vector:"+str(v)+".xyz", disorder=d)
    elif cluster == "FCC":
        for s in sizes:
            i = FCC(atom1="Cu", atom2="Cu", size=s)
            for v in vectors:
                for d in disorder:
                    i.rotate_from_vectors([1, 0, 0], v).get_xyz(file=folder+
                                                                   "/size:"+str(s)+"disorder:"+
                                                                   str(d)+"vector:"+str(v)+".xyz", disorder=d)

    elif cluster == "BCC":
        for s in sizes:
            i = BCC(atom1="Cu", atom2="Cu", size=s)
            for v in vectors:
                for d in disorder:
                    i.rotate_from_vectors([1, 0, 0], v).get_xyz(file=folder+
                                                                   "/size:"+str(s)+"disorder:"+
                                                                   str(d)+"vector:"+str(v)+".xyz", disorder=d)

def simulate_rotate(atom="Cu",
                 folder="all",
                 npt_rot=90,
                 size=1,
                 disorder=0):
    if not os.path.exists(folder):
        os.makedirs(folder)
    cluster = BCC(atom1=atom, atom2=atom, size=size)
    angles = np.linspace(0, np.pi/2, npt_rot)
    vectors =[[np.cos(a), 0, np.sin(a)] for a in angles]
    for v in vectors:
        cluster.rotate_from_vectors([1, 0, 0], v).get_xyz(file=folder +"/vector:" + str(v) + ".xyz", disorder=0.0)


def simulate_series(atom="Cu",
                    folder="all_series",
                    cluster="Icosahedron",
                    sizes=[1, 2, 3, 4],
                    disorder=[0.0, 0.5, 1.0, 1.5, 2.0]):
    if not os.path.exists(folder):
        os.makedirs(folder)
    if cluster == "Icosahedron":
        for s in sizes:
            i = Icosahedron(central_atom="Cu", size=s)
            #5-fold
            #i.rotate_from_vectors(vector1=[0, 0, 1], vector2=[0.0, -0.52573111, 0.85065081], inplace=True)
            #3-fold
            i.rotate_from_vectors(vector1=[0, 0, 1], vector2=[-0.45879397, -0.45879397, 0.45879397], inplace=True)
            for d in disorder:
                i.get_xyz(file=folder+"/size:"+str(s)+"disorder:"+str(d)+".xyz", disorder=d, scale=7)
    elif cluster == "FCC":
        for s in sizes:
            i = FCC(atom1="Cu", atom2="Cu", size=s)
            #3-fold
            i.rotate_from_vectors(vector1=[0, 0, 1], vector2=[1, 1, 0], inplace=True)
            for d in disorder:
                i.get_xyz(file=folder+"/size:"+str(s)+"disorder:"+str(d)+".xyz", disorder=d)
    elif cluster == "BCC":
        for s in sizes:
            i = BCC(atom1="Cu", atom2="Cu", size=s)
            #3-fold
            i.rotate_from_vectors(vector1=[0, 0, 1], vector2=[1, 1, 1], inplace=True)
            for d in disorder:
                i.get_xyz(file=folder+"/size:"+str(s)+"disorder:"+str(d)+".xyz", disorder=d)