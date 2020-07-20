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


def simulate_all_series(folder="all_series",
                        sizes=[3.0, 5.0, 7.0, 9.0, 11.0],
                        disorder=[0.00, 0.20, 0.40, 0.60, 0.80, 1.0],
                        offset=[100, 100, 20]):
    """This function is used to create a series of high symmetry
    atomic positions with varying levels of disorder as well as
    varying size of the cluster
    """
    if not os.path.exists(folder):
         os.makedirs(folder)
    if not os.path.exists(folder + "/BCC6Fold"):
        os.makedirs(folder + "/BCC6Fold")
    if not os.path.exists(folder + "/BCC4Fold"):
        os.makedirs(folder + "/BCC4Fold")
    if not os.path.exists(folder + "/FCC4Fold"):
        os.makedirs(folder + "/FCC4Fold")
    if not os.path.exists(folder + "/FCC6Fold"):
        os.makedirs(folder + "/FCC6Fold")
    if not os.path.exists(folder + "/Ico2Fold"):
        os.makedirs(folder + "/Ico2Fold")
    if not os.path.exists(folder + "/Ico3Fold"):
        os.makedirs(folder + "/Ico3Fold")
    if not os.path.exists(folder + "/Ico5Fold"):
        os.makedirs(folder + "/Ico5Fold")

    for i_s, s in enumerate(sizes):
        # BCC Simulations
        b = BCC(radius=s)
        print("Num BCC Atoms:", len(b))
        b.get_xyz(file=folder + "/BCC4Fold" + "/size:" + str(s), disorder=disorder, offset=offset)
        b.get_6_fold_axis(inplace=True)
        b.get_xyz(file=folder + "/BCC6Fold" + "/size:" + str(s), disorder=disorder, offset=offset)
        # FCC simulations
        f = FCC(radius=s)
        print("Num FCC Atoms:", len(f))
        f.get_xyz(file=folder + "/FCC4Fold" + "/size:" + str(s), disorder=disorder, offset=offset)
        f.get_6_fold_axis()
        f.get_xyz(file=folder + "/FCC6Fold" + "/size:" + str(s), disorder=disorder, offset=offset)
        # Icosahedron Simulations
        i = Icosahedron(num_shells=i_s)
        i.get_xyz(file=folder + "/Ico2Fold" + "/size:" + str(i_s), disorder=disorder, offset=offset)
        i = Icosahedron(num_shells=i_s)
        i.get_3_fold_axis()
        i.get_xyz(file=folder + "/Ico3Fold" + "/size:" + str(i_s), disorder=disorder, offset=offset)
        i.get_5_fold_axis()
        i.get_xyz(file=folder + "/Ico5Fold" + "/size:" + str(i_s), disorder=disorder, offset=offset)
