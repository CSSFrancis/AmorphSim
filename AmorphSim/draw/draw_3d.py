from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
from AmorphSim.utils.vector_utils import get_ico_edges, get_ico_faces


class Icosohedron_3d(Poly3DCollection):
    def __init__(self, line_width=1, alpha=.5, radius=1, location=[0, 0, 0]):
        golden_ratio = .5 * (1 + np.sqrt(5))
        shell_distance = 1
        vertices = [[0, 0, 0],
                    [0, 1, golden_ratio], [0, -1, golden_ratio], [0, 1, -golden_ratio], [0, -1, -golden_ratio],
                    [1, golden_ratio, 0], [-1, golden_ratio, 0], [1, -golden_ratio, 0], [-1, -golden_ratio, 0],
                    [golden_ratio, 0, 1], [-golden_ratio, 0, 1], [golden_ratio, 0, -1], [-golden_ratio, 0, -1]]
        v = [(v / np.linalg.norm(v)) * shell_distance if np.linalg.norm(v) != 0 else v for v in vertices]
        faces, vectors = get_ico_faces(v[1:], shell_dist=shell_distance)
        print("The vectors:", faces)
        new_v = [[v[f[0]],
                  v[f[1]],
                  v[f[2]]] for f in faces]
        new_v = np.multiply(new_v, radius)
        new_v = np.add(new_v, location)
        super().__init__(new_v)
        self.set_alpha(alpha=alpha)
        self.set_linewidth(line_width)
        self.set_edgecolor("k")


class FCC_3d(Poly3DCollection):
    def __init__(self, line_width=1, alpha=.5, radius=1, location= [0, 0, 0]):
        from itertools import permutations
        # create squares
        squares = np.array([[0, 1], [1, 0], [0, -1], [-1, 0]])
        all_squares = []
        for i in range(3):
            temp = squares
            t1 = np.insert(temp, i, values=2, axis=1)
            t2 = np.insert(temp, i, values=-2, axis=1)
            all_squares.append(t1)
            all_squares.append(t2)
        # create hexagons
        base = [[0, -1, 2], [-1, 0, 2], [-2, 0, 1], [-2, -1, 0], [-1, -2, 0], [0, -2, 1]]
        base2 = [[0, -1, -2], [-1, 0, -2], [-2, 0, -1], [-2, -1, 0], [-1, -2, 0], [0, -2, -1]]
        rot90 = [[0, -1, 0], [1, 0, 0], [0, 0, 1]]
        for i in range(4):
            base = np.dot(base, rot90)
            base2 = np.dot(base2, rot90)
            all_squares.append(base)
            all_squares.append(base2)
        all_squares = np.multiply(np.divide(all_squares, np.sqrt(5)),radius)
        print(np.shape(all_squares))
        all_squares = [np.add(i, location) for i in all_squares]
        #new_v = [[v[f[0]], v[f[1]], v[f[2]]] for f in faces]
        super().__init__(all_squares)
        self.set_alpha(alpha=alpha)
        self.set_linewidth(line_width)
        self.set_edgecolor("k")


class Cube3d(Poly3DCollection):
    def __init__(self, line_width=1, alpha=.5, shape= [20,20,20], location= [0, 0, 0]):
        # create squares
        v = np.array([[0, 0, 0],
                      [shape[0], 0,0],
                      [shape[0], shape[1], 0],
                      [0, shape[1], 0],
                      [shape[0], 0, shape[2]],
                      [0, 0, shape[2]],
                      [0, shape[1], shape[2]],
                      [shape[0],shape[1], shape[2]]])

        verticies = [[v[0], v[1], v[2], v[3]],
                     [v[4], v[5], v[6], v[7]],
                     [v[0], v[1], v[4], v[5]],
                     [v[2], v[3], v[6], v[7]],
                     [v[2], v[1], v[4], v[7]],
                     [v[5], v[6], v[3], v[0]]]

        super().__init__(verticies)
        self.set_alpha(alpha=alpha)
        self.set_linewidth(line_width)
        self.set_edgecolor("k")
        self.set_facecolor("k")


