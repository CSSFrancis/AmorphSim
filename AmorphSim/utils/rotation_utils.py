import numpy as np
from numpy.random import random

def _get_random_rotation(num_vectors):
    rand_vector = random((num_vectors, 3)) * 2 - 1  # change
    norms = np.linalg.norm(rand_vector, axis=1)
    rand_vector = np.divide(rand_vector,norms[:,np.newaxis])
    rand_rot = random(num_vectors) * np.pi
    return rand_vector, rand_rot

def _rand_2d_rotation_matrix():
    angle = np.random.uniform(low=0.0, high=np.pi)
    M = np.array([[np.cos(angle),-np.sin(angle), 0],
                  [np.sin(angle), np.cos(angle), 0],
                  [0,0,1]])
    return M


def _rand_3d_rotation_matrix(deflection=1.0, randnums=None):
    """
    Creates a random rotation matrix.

    deflection: the magnitude of the rotation. For 0, no rotation; for 1, competely random
    rotation. Small deflection => small perturbation.
    randnums: 3 random numbers in the range [0, 1]. If `None`, they will be auto-generated.
    """
    # from http://www.realtimerendering.com/resources/GraphicsGems/gemsiii/rand_rotation.c
    randnums = np.random.uniform(size=(3,))
    theta, phi, z = randnums
    theta = theta * 2.0 * deflection * np.pi  # Rotation about the pole (Z).
    phi = phi * 2.0 * np.pi  # For direction of pole deflection.
    z = z * 2.0 * deflection  # For magnitude of pole deflection.
    r = np.sqrt(z)
    Vx, Vy, Vz = V = (np.sin(phi) * r, np.cos(phi) * r, np.sqrt(2.0 - z))
    st = np.sin(theta)
    ct = np.cos(theta)
    R = np.array(((ct, st, 0), (-st, ct, 0), (0, 0, 1)))
    # Construct the rotation matrix  ( V Transpose(V) - I ) R.
    M = (np.outer(V, V) - np.eye(3)).dot(R)
    return M


def _get_random_2d_rot(num_vectors):
    rot_m = []
    for i in range(0, num_vectors):
        M = _rand_2d_rotation_matrix()
        rot_m.append(M)
    return rot_m


def _get_random_3d_rot(num_vectors):
    rot_m = []
    for i in range(0, num_vectors):
        M = _rand_3d_rotation_matrix()
        rot_m.append(M)
    return rot_m


def _get_rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians. Using Euler-Rodrigues Formula
    """
    axis = np.asarray(axis)
    axis = axis / np.sqrt(np.dot(axis, axis))
    a = np.cos(theta / 2.0)
    b, c, d = -axis * np.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])



