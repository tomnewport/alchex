from cgswap.geometry import PointCloud, plot_3d

import numpy

a = PointCloud(3)
a.add_points(numpy.random.randn(10,3))

b = a.interpolate_1d_list(numpy.linspace(0,1,100))

plot_3d(a, b)

