import numpy as np
#import hyperspy.api as hs


def _get_speckle_intensity(k_vector,
                           ewald_sphere_rad,
                           cluster_rad = 1.0,
                           disorder=None,
                           beam_direction=[0,0,1]):
    """Returns the intensity of some speckle out of one

    Parameters
    --------------
    k_vector:
        The vector for one of the diffraction spots
    ewald_sphere_rad:
        The radius of the Ewald Sphere.  This is equal to 1/wavelength
    disorder: None or float
        The amount of disorder in the glass.  This is a measure of the average
        displacement of some atom from ideality.
    :returns
    -----------
    observed_intensity:
        The intensity for some diffraction vector
    :return:
    """
    deviation = _get_deviation(ewald_sphere_rad,
                               k_vector,
                               beam_direction=beam_direction)
    observed_intensity =_shape_function(radius=cluster_rad,
                                        deviation=deviation)
    if disorder is not None:
        factor = _get_disorder(k_vector,disorder)
        observed_intensity = observed_intensity*factor
    return observed_intensity

def _shape_function(radius, deviation, function=None):
    """Returns the point at some deviation parameter s and a radius r

    Parameters
    ----------------
    r: float
        The radius of the cluster
    s_g: float
        The deviation parameter for the point.
    function: None
        The shape function to simulate (still need to add different functions)

    Returns:
    n:float
        Shape function for some deviation parameter and radius
    """
    c = 2* np.pi*deviation*radius
    n = ((3*(np.sin(c)-c*np.cos(c)))/c**3)**2
    return n


def s_g_kernel(kernel_size, d_hkl, cluster_size, voltage):
    """ Simulates a 2-d projection of the s_g kernel... (Maybe make 3-D if useful)
    Parameters
    ----------------
    kernel_size: int
        The size of the kernel to create
    d_hkl: float
        The interplanar spacing in n,
    cluster_size: float
        The size of the cluster being calculated.
    """
    wavelength = _get_wavelength(acc_voltage=voltage)
    ax = np.arange(-kernel_size // 2 + 1., kernel_size // 2 + 1.)
    xx, yy = np.divide(np.meshgrid(ax, ax), kernel_size / 3 * cluster_size)
    scaling = 1/(kernel_size / 3 * cluster_size)
    sg = np.power((np.square(xx) + np.square(yy)), 0.5)
    print(d_hkl)
    dot = np.subtract(2*wavelength*d_hkl, np.multiply(sg, d_hkl))
    angles = np.divide(sg, d_hkl)
    sg_surf = np.multiply(sg, (2 * np.pi * cluster_size))
    kernel = np.power(np.multiply(np.divide((np.sin(sg_surf) -
                                             np.multiply(sg_surf,
                                                         np.cos(sg_surf))),
                                            np.power(sg_surf, 3)), 3), 2)
    dict0 = {'size': kernel_size, 'name': 's_x', 'units': 'nm^-1', 'scale': scaling, 'offset': 0}
    dict1 = {'size': kernel_size, 'name': 's_y', 'units': 'nm^-1', 'scale': scaling, 'offset': 0}
    k = hs.signals.Signal2D(data=kernel, axes=[dict0, dict1])
    return k

def s_g_kernel_angle(kernel_size, d_hkl, cluster_size, voltage):
    """ Simulates a 2-d projection of the s_g kernel... (Maybe make 3-D if useful)
    Parameters
    ----------------
    kernel_size: int
        The size of the kernel to create
    d_hkl: float
        The interplanar spacing in n,
    cluster_size: float
        The size of the cluster being calculated.
    """
    wavelength = _get_wavelength(acc_voltage=voltage)
    ax = np.arange(-kernel_size // 2 + 1., kernel_size // 2 + 1.)
    xx, yy = np.divide(np.meshgrid(ax, ax), kernel_size / 3 * cluster_size)
    scaling = 1/(kernel_size / 3 * cluster_size)
    sg = np.power((np.square(xx) + np.square(yy)), 0.5)
    print(d_hkl)
    dot = np.subtract(2*wavelength*d_hkl, np.multiply(sg, d_hkl))
    angles = np.divide(sg, d_hkl)
    sg_surf = np.multiply(sg, (2 * np.pi * cluster_size))
    kernel = np.power(np.multiply(np.divide((np.sin(sg_surf) -
                                             np.multiply(sg_surf,
                                                         np.cos(sg_surf))),
                                            np.power(sg_surf, 3)), 3), 2)
    dict0 = {'size': kernel_size, 'name': 's_x', 'units': 'nm^-1', 'scale': scaling, 'offset': 0}
    dict1 = {'size': kernel_size, 'name': 's_y', 'units': 'nm^-1', 'scale': scaling, 'offset': 0}
    k = hs.signals.Signal2D(data=kernel, axes=[dict0, dict1])
    return k


def _get_speckle_size(accelerating_voltage=200, semi_angle= 0.74):
    """Gets the size of a speckle from the semi angle
    Parameters
    ----------------
    accelerating_voltage:
        The accelerating voltage for the microscope
    semi_angle: float
        The angle in mili-radians of the convergance
    k: float
        The k-vector length in nm
    """
    wavelength = _get_wavelength(acc_voltage=accelerating_voltage)
    size = np.sin(semi_angle/1000) * 1/wavelength
    return size


def _get_wavelength(acc_voltage):
    """Given some accelerating voltage for a microscope calculate the relativistic wavelength
    Parameters
    -----------------
    acc_voltage: float
        The accelerating voltage of the microscope.

    Returns:
    -----------
    wavelength:float
        The wavelength of the electrons.
    """
    h = 6.626*10**-34
    m0 = 9.109*10**-31
    e = 1.602*10**-19
    wavelength = h/np.sqrt(2*m0*acc_voltage*1000*e*(1+acc_voltage/(2*511)))*10**9
    return wavelength


def _get_deviation(sphere_radius, k, beam_direction=[0,0,-1]):
    """
    Parameters
    ----------------
    sphere_radius: float
        The radius of the sphere
    k0: tuple
        The (x,y,z) of the original s value from the optic axis.
    """
    beam = np.array(beam_direction) * sphere_radius
    dist = np.linalg.norm(beam-k)
    deviation = sphere_radius-dist # distance from edge of sphere to k
    return deviation

def _get_disorder(k, displacement):
    """The estimates the amount of disorder in the perfect crystal.
    This might need some additional work.

    Parameters
    ----------------
    sphere_radius: float
        The radius of the sphere
    k: tuple
        The (x,y,z) of the original s value from the optic axis.
    """
    b = 8 * np.pi**2 * displacement**2
    s = np.linalg.norm(k)
    factor = np.exp(-b*s**2)
    return factor