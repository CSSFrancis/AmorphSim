import numpy as np


def _get_deviation(sphere_radius,k):
    """
    Parameters
    ----------------
    sphere_radius: float
        The radius of the sphere
    k0: tuple
        The (x,y,z) of the original s value from the optic axis.
    """
    dist = np.sqrt(k[0]**2+k[1]**2+(-sphere_radius - k[2])**2) # Distance from the center of sphere to k
    deviation = sphere_radius-dist # distance from edge of sphere to k
    return deviation


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



