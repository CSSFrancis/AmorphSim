
class GlassCluster(object):
    def __init__(self, vertexes=None, radius=1, k=4.0, position=random(2),
                 rotation_vector=[1, 0, 0], rotation_angle=0, diffraction_intensity=600):
        """Defines a cluster with a symmetry of symmetry, a radius of radius in nm and position of position.

        Parameters:
        ----------------
        symmetry: int
            The symmetry of the cluster being simulated
        radius: float
            The radius of the cluster in nm
        position: tuple
            The position of the cluster in the simulation cube
        rotation_vector: tuple
            The vector which the cluster is rotated around
        rotation_angle: float
            The angle the cluster is rotated about
        """
        self.symmetry = symmetry
        self.radius = radius
        self.position = position
        self.rotation_vector = rotation_vector
        self.rotation_angle = rotation_angle
        self.k = k
        self.diffraction_intensity = diffraction_intensity

    def get_diffraction(self, img_size=8.0, num_pixels=512, accelerating_voltage=200, conv_angle=0.6):
        """Takes some image size in inverse nm and then plots the resulting
        """
        rotation_matrix = _get_rotation_matrix(self.rotation_vector, self.rotation_angle)
        sphere_radius = 1/_get_wavelength(accelerating_voltage)
        scale = (num_pixels-1)/img_size
        angle = (2 * np.pi) / self.symmetry  # angle between speckles on the pattern
        k = [[np.cos(angle * i) * self.k, np.sin(angle * i) * self.k,0] for i in
             range(self.symmetry)]  # vectors for the speckles perp to BA
        k_rotated = [np.dot(rotation_matrix, speckle) for speckle in k]
        deviation = [_get_deviation(sphere_radius,speckle) for speckle in k_rotated]
        observed_intensity = [self.diffraction_intensity/self.symmetry * _shape_function(radius=self.radius, deviation=dev)
                              for dev in deviation]
        radius = _get_speckle_size(accelerating_voltage, conv_angle)*scale
        circles = [circle(int(k1[0] * scale + num_pixels/2), int(k1[1] * scale + num_pixels/2),
                          radius=radius) for k1 in k_rotated]
        image = np.ones(shape=(num_pixels,num_pixels))
        for (r, c), i in zip(circles, observed_intensity):
            image[r, c] = i + image[r, c]
        image = gaussian(image=image, sigma=2)
        return image

    def get_speckles(self, img_size=8.0, num_pixels=512, accelerating_voltage=200, conv_angle=0.6):
        """Takes some image size in inverse nm and then plots the resulting
        """
        rotation_matrix = _get_rotation_matrix(self.rotation_vector, self.rotation_angle)
        sphere_radius = 1/_get_wavelength(accelerating_voltage)
        scale = (num_pixels-1)/img_size
        angle = (2 * np.pi) / self.symmetry  # angle between speckles on the pattern
        k = [[np.cos(angle * i) * self.k, np.sin(angle * i) * self.k,0] for i in
             range(self.symmetry)]  # vectors for the speckles perp to BA
        k_rotated = [np.dot(rotation_matrix, speckle) for speckle in k]
        deviation = [_get_deviation(sphere_radius,speckle) for speckle in k_rotated]
        observed_intensity = [self.diffraction_intensity/self.symmetry * _shape_function(radius=self.radius, deviation=dev)
                              for dev in deviation]
        radius = _get_speckle_size(accelerating_voltage, conv_angle)*scale
        speckles = [circle(int(k1[0] * scale + num_pixels/2), int(k1[1] * scale + num_pixels/2),
                           radius=radius) for k1 in k_rotated]
        return speckles, observed_intensity

    def get_intensity(self,accelerating_voltage=200):
        """Takes some image size in inverse nm and then plots the resulting
        """
        rotation_matrix = _get_rotation_matrix(self.rotation_vector, self.rotation_angle)
        sphere_radius = 1/_get_wavelength(accelerating_voltage)
        angle = (2 * np.pi) / self.symmetry  # angle between speckles on the pattern
        k = [[np.cos(angle * i) * self.k, np.sin(angle * i) * self.k,0] for i in
             range(self.symmetry)]  # vectors for the speckles perp to BA
        k_rotated = [np.dot(rotation_matrix, speckle) for speckle in k]
        deviation = [_get_deviation(sphere_radius,speckle) for speckle in k_rotated]
        observed_intensity = [self.diffraction_intensity/self.symmetry * _shape_function(radius=self.radius, deviation=dev)
                              for dev in deviation]
        return observed_intensity

    def get_mistilt(self, v2=[1,1,1]):
        v1 = self.rotation_vector
        print(v1)
        return np.cos((v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])/(np.linalg.norm(v1)*np.linalg.norm(v2)))