from AmorphSim.create_xyz import simulate_all, simulate_rotate, simulate_series
from unittest import TestCase



class TestXYZ(TestCase):
    def test_simualte_all(self):
        simulate_all(cluster="BCC", folder="AllBCC")

    def test_simualte_roations(self):
        simulate_rotate(size=1,folder="RotSize1BCC")

    def test_simulate_series(self):
        simulate_series(cluster="Icosahedron",folder="IcoSeries3fold")

    def test_simulate_series2(self):
        simulate_series(cluster="BCC", folder="BCCSeries6fold")

    def test_simulate_series3(self):
        simulate_series(cluster="FCC", folder="FCCSeries4fold")