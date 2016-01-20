import numpy
from colorsys import hsv_to_rgb
import matplotlib.pyplot as plt
from random import randint, random
from mpl_toolkits.mplot3d import Axes3D

class TransformationMatrix:
    def __init__(self, dimensions=2):
        self.dimensions = dimensions
        self.matrix = numpy.zeros((dimensions+1, dimensions+1))
        for x in range(dimensions+1):
            self.matrix[x,x] = 1
    def translate(self, *axes):
        sxes = axes[:self.dimensions]
        transmatrix = numpy.zeros((self.dimensions+1, self.dimensions+1))
        numpy.fill_diagonal(transmatrix, 1)
        for axis, axistranslation in enumerate(axes):
            transmatrix[axis,-1] += axistranslation
        self.matrix = self.matrix.dot(transmatrix)
    def randomise_3d(self):
        self.rotate_3d(random(), [random(), random(), random()])
        self.translate(random(), random(), random())
    def rotate_3d(self, theta, vector):
        vector = numpy.array(vector)/numpy.linalg.norm(vector)
        l = vector[0]
        m = vector[1]
        n = vector[2]
        cos_theta = numpy.cos(theta)
        sin_theta = numpy.sin(theta)
        self.matrix = self.matrix.dot(numpy.array([
            [
                l*l*(1-cos_theta) + cos_theta, 
                m*l*(1-cos_theta) - (n * sin_theta),
                n*l*(1-cos_theta) + (m * sin_theta),
                0
            ],
            [
                l*m*(1-cos_theta) + n * sin_theta,
                m*m*(1-cos_theta) + cos_theta,
                n*m*(1-cos_theta) - l*sin_theta,
                0
            ],
            [
                l*n*(1-cos_theta) - m*sin_theta,
                m*n*(1-cos_theta) + l*sin_theta,
                n*n*(1-cos_theta) + cos_theta,
                0
            ],
            [
                0,
                0,
                0,
                1
            ]
        ]))
    def rotate_3d_about(self, point, theta, axis):
        self.translate(*point)
        self.rotate_3d(theta, axis)
        self.translate(*[-x for x in point])
    def rotate_2d(self, theta, vector):
        x = vector[0]
        y = vector[1]
        cos_theta = numpy.cos(theta)
        sin_theta = numpy.sin(theta)
        self.matrix = self.matrix.dot(numpy.array(
            [[cos_theta, -sin_theta, x - (cos_theta*x) - (-sin_theta*y)],
             [sin_theta, cos_theta, y-(sin_theta*x) - (cos_theta*y)],
             [0,0,1]]).astype(numpy.float64))

class PointCloud:
    def __init__(self, dimensions=2):
        self.dimensions = dimensions
        self.points = numpy.zeros((0,self.dimensions+1))
        self.lines  = []
    def randomsetup(self, points=10, lines=10, size=1):
        self.add_points(numpy.random.randn(points,self.dimensions) * size)
        attempts = 0
        while len(self.lines) < lines and attempts < lines * 10:
            attempts += 1
            line = [randint(0,points-1), randint(0,points-1)]
            if line not in self.lines and line[::-1] not in self.lines and len(set(line)) == 2:
                self.lines.append(line)
    def add_points(self, coordinates):
        coordinates = numpy.pad(
            coordinates,
            pad_width= ((0,0),(0,1)), 
            mode="constant", 
            constant_values=1
        )
        self.points = numpy.append(self.points, coordinates, axis=0)
    def add_line(self, point1, point2):
        self.lines.append((point1, point2))
    def mutate(self, size=0.1):
        self.points[:,:3] += (numpy.random.randn(*self.points.shape) * size)[:,:3]
    def _plot_2d(self):
        for line_id, (point1, point2) in enumerate(self.lines):
            h = (line_id*0.5)/len(self.lines)
            coordinates_1 = self.points[point1,:]
            coordinates_2 = self.points[point2,:]
            plt.plot(
                [coordinates_1[0],coordinates_2[0]], 
                [coordinates_1[1],coordinates_2[1]], 
                color=hsv_to_rgb(h,1,1))
        plt.scatter(self.points[:,0],self.points[:,1])
        plt.axis('equal')
    def _plot_3d(self, *transformation_matrices):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        hues = numpy.linspace(0,0.8,self.points.shape[0])
        transformation_matrices = list(transformation_matrices)
        transformation_matrices.append(None)
        plottable = self.clone()
        for transformation_matrix, l in zip(transformation_matrices, numpy.linspace(1,0.1,len(transformation_matrices))):
            cols = [hsv_to_rgb(h,l,1) for h in hues]
            if transformation_matrix is not None:
                plottable.transform(transformation_matrix)
            ax.scatter(
                plottable.points[:,0], 
                plottable.points[:,1], 
                plottable.points[:,2], 
                 lw = 0,
                alpha=0.7,
                c=cols,
                s=50
            )
    def interpolate(self, other, l=0.5):
        interpolated = self.clone()
        interpolated.points = other.points * l + self.points * (1-l)
        return interpolated
    def plot_3d_transform(self, transformation, references=(), steps = 10):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        hues = numpy.linspace(0,0.8,self.points.shape[0])
        transformed = self.clone().transform(transformation)
        if len(references) > 0:
            a = 0.2
        else:
            a = 0.5
        for l in numpy.linspace(0,1,steps):
            interpolated = self.interpolate(transformed, l)
            cols = [hsv_to_rgb(h,(0.5+l/2),1) for h in hues]
            lw = 1 if l in [1,0] else 0
            ax.scatter(
                interpolated.points[:,0], 
                interpolated.points[:,1], 
                interpolated.points[:,2], 
                lw = lw,
                alpha=a,
                c=cols,
                s=50
            )
        for r in references:
            cols = [hsv_to_rgb(h,1,1) for h in hues]
            lw = 1 if l in [1,0] else 0
            ax.scatter(
                r.points[:,0], 
                r.points[:,1], 
                r.points[:,2], 
                lw = 1,
                alpha=0.8,
                c=cols,
                s=100
            )
    def paired_3d_align(self, other, inv=True, centroid_weighting=None):
        # inv = True  - move this pointcloud to the other
        #       False - move the other pointcloud to this
        if inv:
            return solve_3d_transformation(self, other, centroid_weighting=centroid_weighting)
        else:
            return solve_3d_transformation(other, self, centroid_weighting=centroid_weighting)
    def n_points(self):
        return self.points.shape[0]
    def unpaired_3d_rmse(self, other):
        if self.n_points < other.n_points:
            shorter = self
            longer = other
        else:
            longer = self
            shorter = other
        return distance.cdist(longer, shorter).min(axis=0).mean()
    def centroid(self):
        return numpy.mean(self.points, axis=0)
    def bounding_box(self, extent=True):
        mins = self.points.min(axis=0)
        maxs = self.points.max(axis=0)
        if extent:
            centroid = self.centroid()
            return mins-centroid, maxs-centroid
        else:
            return mins, maxs
    def plot(self, *args):
        if self.dimensions == 2:
            self._plot_2d()
        if self.dimensions == 3:
            self._plot_3d(*args)
    def transform(self, transformation_matrix):
        self.points = self.points.astype(numpy.float64)
        self.points = self.points.dot(transformation_matrix.matrix.T)
        return self
    def clone(self):
        clone = PointCloud(self.dimensions)
        clone.points = self.points
        clone.lines = self.lines
        return clone
    
def pairwise_rmse(points1, points2):
    return numpy.sqrt((numpy.square(points1 - points2)).mean())

def solve_3d_transformation(operand, reference, subselection=None, test=True, centroid_weighting=None):
    if subselection is None:
        operand_points = numpy.mat(operand.points[:,:3])
        reference_points = numpy.mat(reference.points[:,:3])
    else:
        operand_points = numpy.mat(operand.points[subselection,:3])
        reference_points = numpy.mat(reference.points[subselection,:3])
    assert len(operand_points) == len(reference_points)

    n_points = operand_points.shape[0]
    print(centroid_weighting)
    if centroid_weighting is None:
        centroid_weighting = [1] * n_points
    operand_centroid = numpy.average(operand_points, axis=0, weights=centroid_weighting)
    reference_centroid = numpy.average(reference_points, axis=0, weights=centroid_weighting)
    
    # centre the points
    operand_centred   = operand_points - numpy.tile(operand_centroid, (n_points, 1))
    reference_centred   = reference_points - numpy.tile(reference_centroid, (n_points, 1))

    # dot is matrix multiplication for array
    covariance_matrix = numpy.transpose(operand_centred) * reference_centred

    unitary_matrix, diagonal_matrix, unitary_over_field = numpy.linalg.svd(covariance_matrix)

    rotation_matrix = unitary_over_field.T * unitary_matrix.T
    # special reflection case
    if numpy.linalg.det(rotation_matrix) < 0:
       unitary_over_field[2,:] *= -1
       rotation_matrix = unitary_over_field.T * unitary_matrix.T
    translation_vector = -rotation_matrix*operand_centroid.T + reference_centroid.T

    rotation_matrix = numpy.pad(
            rotation_matrix, 
            pad_width = ((0,1),(0,1)), 
            mode="constant", 
            constant_values=0)
    rotation_matrix[3,3] = 1


    translation_matrix = numpy.zeros((4,4))
    numpy.fill_diagonal(translation_matrix, 1)
    translation_matrix[:3,3:] = translation_vector
    

    transformation_matrix = TransformationMatrix(3)
    transformation_matrix.matrix = translation_matrix.dot(rotation_matrix)

    if test:
        operand_4d =numpy.pad(
            operand_points,
            pad_width= ((0,0),) * (1) + ((0,1),), 
            mode="constant", 
            constant_values=1
        )
        #transformed_operand = operand_4d.dot(transformation_matrix.matrix.T)[:,:3]
        transformed_operand = operand.clone()
        transformed_operand.transform(transformation_matrix)
        rmse = pairwise_rmse(transformed_operand.points[:,:3], reference.points[:,:3])
        #print(rmse)
        return transformation_matrix, rmse, transformed_operand
    else:
        return transformation_matrix


def plot_3d(*pointclouds):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_aspect('equal')
    hues = numpy.linspace(0,0.8,len(pointclouds))
    xmin = 10000000
    ymin = 10000000
    zmin = 10000000
    xmax = -10000000
    ymax = -10000000
    zmax = -10000000
    for hue, plottable in zip(hues, pointclouds):
        cols = hsv_to_rgb(hue,1,1)
        ax.scatter(
                plottable.points[:,0], 
                plottable.points[:,1], 
                plottable.points[:,2], 
                 lw = 0,
                alpha=0.7,
                c=cols,
                s=100
            )
        xmin = min(xmin, plottable.points[:,0].min())
        xmax = max(xmax, plottable.points[:,0].max())
        ymin = min(ymin, plottable.points[:,1].min())
        ymax = max(ymax, plottable.points[:,1].max())
        zmin = min(zmin, plottable.points[:,2].min())
        zmax = max(zmax, plottable.points[:,2].max())
    box = max(xmax - xmin, ymax - ymin, zmax - zmin)
    xmid = (xmax + xmin) / 2.0
    ymid = (ymax + ymin) / 2.0
    zmid = (zmax + zmin) / 2.0
    ax.scatter(
        [xmid-box/1.8, xmid+box/1.8, xmid-box/1.8, xmid+box/1.8, xmid-box/1.8, xmid+box/1.8, xmid-box/1.8, xmid+box/1.8],
        [ymid-box/1.8, ymid-box/1.8, ymid+box/1.8, ymid+box/1.8, ymid-box/1.8, ymid-box/1.8, ymid+box/1.8, ymid+box/1.8],
        [zmid-box/1.8, zmid-box/1.8, zmid-box/1.8, zmid-box/1.8, zmid+box/1.8, zmid+box/1.8, zmid+box/1.8, zmid+box/1.8])
    ax.set_aspect('equal')
    plt.show()