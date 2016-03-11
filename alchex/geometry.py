import numpy
from colorsys import hsv_to_rgb
import matplotlib.pyplot as plt
from random import randint, random
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial.distance import pdist, cdist, squareform
#import seaborn as sns
import networkx as nx
import matplotlib.pyplot as plt
from math import floor
from alchex.algorithms import bruteforce_pair

def linear_weights(coordinate, length):
    w = numpy.zeros(length)
    v_index = coordinate * (length-1)
    f = int(floor(v_index))
    c = f+1
    w[f] = c - v_index
    if c <= (length-1):
        w[c] = v_index - f
    return w

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
    def scale(self, *axes):
        scalematrix = numpy.zeros([self.dimensions+1, self.dimensions+1])
        axes = list(axes) + [1]
        numpy.fill_diagonal(scalematrix, axes)
        self.matrix = self.matrix.dot(scalematrix)
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
        self.faces = []
        self.vertex_normals = None
    def group_points(self, group_size, **kwargs):
        if group_size == 2:
            return bruteforce_pair(self.points, **kwargs)
        else:
            raise NotImplementedError()
    def find_friends(self, other, n_friends, distance_cutoff=15):
        distance_matrix = cdist(self.points, other.points)
        distance_matrix[distance_matrix==0] = 1000000000
        distance_ranks = numpy.argsort(-distance_matrix, axis=None)
        #sns.distplot(distance_matrix.min(axis=0))
        #sns.plt.show()
        friends = []
        used = set()
        depth = 0
        while len(friends) < n_friends and depth < distance_matrix.shape[0] ** 2:
            depth += 1
            p1, p2 = numpy.unravel_index(distance_ranks[-depth], distance_matrix.shape)
            #p1, p2 = pair[0]
            if p1 not in used and p2 not in used:
                used.update({p1, p2})
                friends.append((p1, p2, distance_matrix[p1, p2]))
        return friends
    def interpolate_1d(self, coordinate, weighting_function="linear"):
        if weighting_function == "linear":
            weights = linear_weights(coordinate, self.points.shape[0])
        return self.centroid(centroid_weighting=weights)
    def interpolate_1d_list(self, coordinates, weighting_function="linear"):
        r = PointCloud(self.dimensions)
        r.points = numpy.array([self.interpolate_1d(c, weighting_function=weighting_function) for c in coordinates])
        return r
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
        for line in self.lines:
            p1 = self.points[line[0],:]
            p2 = self.points[line[1],:]
            ax.plot([p1[0],p2[0]], [p1[1],p2[1]], [p1[2],p2[2]])
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
    def centroid(self, centroid_weighting=None):
        if centroid_weighting is None:
            centroid_weighting = [1]*self.points.shape[0]
        return numpy.average(self.points, axis=0, weights=centroid_weighting)
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
    def cliques_to_faces(self, max_faces=3):
        graph = nx.Graph()
        graph.add_edges_from(self.lines)
        print("finding")
        cliques = sorted(nx.find_cliques(graph), key=lambda x: -len(x))
        self.faces = []
        for clique in cliques:
            h_clique = frozenset(clique)
            if len(clique) == 3:
                self.faces.append(clique)
            else:
                clique_positions = self.points[clique,:3]
                ext_edges = []
                int_edges = []
                for edge in combinations(clique,2):
                    shared_triangles = set([x[1] for x in graph.edges(edge[0])]).intersection(set([x[1] for x in graph.edges(edge[1])]))
                    print shared_triangles
                    ext_connections = shared_triangles - set(clique)
                    print ext_connections
                    if len(ext_connections) > 0:
                        ext_edges.append(edge)
                    else:
                        int_edges.append(edge)
                print len(ext_edges), len(int_edges)

class Volume(PointCloud):
    pass




        
    def build_triangles(self):
        '''
        This is a fast way to skin the mesh assuming
        there are no cliques of size > 3.
        '''
        graph = nx.Graph()
        graph.add_edges_from(self.lines)
        for n1 in graph:
            for e2 in graph.edges(n1):
                n2 = e2[1]
                if n2 == n1:
                    continue
                for e3 in graph.edges(n2):
                    n3 = e3[1]
                    if n3 == n2:
                        continue
                    for e4 in graph.edges(n3):
                        n4 = e4[1]
                        if n4 == n1:
                            self.faces.append([n1, n2, n3])

    
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

from scipy.spatial import Delaunay, ConvexHull
from itertools import combinations, product
from collections import Counter

def export_obj_3d(pointcloud, destination):
    with open(destination, "w") as obj_handle:
        for point in pointcloud.points:
            obj_handle.write("v "+ " ".join(map(str,point))+"\n")
        if pointcloud.vertex_normals is not None:
            for vert_norm in pointcloud.vertex_normals:
                obj_handle.write("vn "+" ".join(map(str, vert_norm))+"\n")
        for line in pointcloud.lines:
            obj_handle.write("l "+ " ".join(map(lambda x: str(x+1),line))+"\n")
        for plane in pointcloud.faces:
            obj_handle.write("f " + " ".join(map(lambda x : str(x+1), plane))+"\n")


def alpha_shape_3d(pointcloud, alpha=1):
    all_simplices = Delaunay(pointcloud.points[:,:3]).simplices
    matching_simplices = []
    faces = []
    all_points = set(range(pointcloud.points.shape[0]))
    simplex_points = set()
    for simp in all_simplices:
        simp_radius = 0
        for a, b in combinations(simp,2):
            simp_radius = max(simp_radius, numpy.linalg.norm(pointcloud.points[a,:] - pointcloud.points[b,:]))
        if simp_radius < alpha:
            simplex_points.update(simp)
            matching_simplices.append(simp)
    for simp in matching_simplices:
        faces += [frozenset(x) for x in combinations(simp,3)]
    outer = [plane for plane, count in Counter(faces).items() if count == 1]
    hull_points = [x for y in outer for x in y]
    for plane in outer:
        for a, b in combinations(plane,2):
            pointcloud.lines.append([a,b])
        pointcloud.faces.append(list(plane))
    lone_points = all_points - simplex_points
    export_obj_3d(pointcloud, "test.obj")

def voronoi_shell_3d(pointcloud, subpoints):
    subpoints = set(subpoints)
    print("dl-doing")
    all_simplices = Delaunay(pointcloud.points[:,:3]).simplices
    print("dl-done")
    shell_simplices = []
    shell_nodes = set()
    shell_node_norms = {}
    shell_vnorms = {}
    shell_edges = set()
    shell_faces = []
    for simplex in all_simplices:
        simplex = set(simplex)
        simplex_subnodes = simplex.intersection(subpoints)
        simplex_supnodes = simplex - simplex_subnodes
        simplex_edges = [x for x in product(simplex_subnodes,simplex_supnodes)]
        simplex_faces = map(frozenset,(combinations(simplex, 3)))
        if len(simplex_edges) >= 3:
            shell_nodes.update(simplex_edges)
            for inside_node, outside_node in simplex_edges:
                if (inside_node, outside_node) not in shell_node_norms:
                    # Normal is pointing outside, so inside -> outside is 
                    # the right vector
                    shell_node_norms[(inside_node, outside_node)] = pointcloud.points[outside_node,:3] - pointcloud.points[inside_node,:3]
                    shell_vnorms[(inside_node, outside_node)] = []
            this_shell_edges = {x:[] for x in simplex_edges}
            for e1, e2 in combinations(simplex_edges,2):
                i = set(e1).intersection(set(e2))
                if len(i) == 1:
                    this_shell_edges[e1].append(e2)
                    this_shell_edges[e2].append(e1)
                    shell_edges.add((e1,e2))
            this_shell_face = [this_shell_edges.keys()[0]]
            while len(this_shell_face) < len(this_shell_edges):
                next_nodes = this_shell_edges[this_shell_face[-1]]
                for node in next_nodes:
                    if node not in this_shell_face:
                        this_shell_face.append(node)
                        break
            face_points = [numpy.mean(pointcloud.points[edge,:3], axis=0) for edge in this_shell_face]
            face_norms  = [shell_node_norms[x] for edge in this_shell_face]
            face_norm = numpy.mean(face_norms, axis=0)
            pvec1 = face_points[0] - face_points[1]
            pvec2 = face_points[1] - face_points[2]
            plane_normal = numpy.cross(pvec1, pvec2)
            vec_normal = numpy.mean(face_norms, axis=0)
            face_centre = numpy.mean(face_points, axis=0)
            # Check and flip
            if numpy.dot(plane_normal, vec_normal) < 0:
                plane_normal = - plane_normal
                this_shell_face = this_shell_face[::-1]
            shell_faces.append(this_shell_face)
            for edge in this_shell_face:
                shell_vnorms[edge].append(plane_normal)
    shell_nodes = list(shell_nodes)
    shell_lookup = {x:idx for idx, x in enumerate(shell_nodes)}
    shell_vnorms = [shell_vnorms[x] for x in shell_nodes]
    shell_nodes = [pointcloud.points[list(x),:3].mean(axis=0) for x in shell_nodes]
    a = Volume(3)
    a.add_points(shell_nodes)
    #a.vertex_normals = shell_node_norms
    for n1, n2 in shell_edges:
        a.lines.append([shell_lookup[n1], shell_lookup[n2]])
    for face in shell_faces:
        a.faces.append([shell_lookup[x] for x in face])
    export_obj_3d(a, "test.obj")





def cross_sectional_area_3d(pointcloud, axis, method="convexhull"):
    plane_normal = numpy.cross(axis, [0,0,1])
    rot_angle = numpy.arctan2(numpy.linalg.norm(plane_normal), numpy.dot(axis,[0,0,1]))

    projection = TransformationMatrix(3)
    projection.rotate_3d(rot_angle, plane_normal)

    projected = pointcloud.clone()
    projected.transform(projection)
    projected.points[:,2] = 1
    projected.dimensions = 2
    projected.points = projected.points[:,:3]
    if pointcloud.points.shape[0] < 3:
        return 0
    if method == "convexhull":
        triangles = Delaunay(projected.points[:,:2]).simplices
        area = 0
        for tp1, tp2, tp3 in triangles:
            tl12 = numpy.linalg.norm(projected.points[tp1,:2] - projected.points[tp2,:2])
            tl23 = numpy.linalg.norm(projected.points[tp2,:2] - projected.points[tp3,:2])
            tl34 = numpy.linalg.norm(projected.points[tp3,:2] - projected.points[tp1,:2])
            s = (tl12 + tl23 + tl34) / 2
            area += (s*(s-tl12)*(s-tl23)*(s-tl34)) ** 0.5
        return area

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
                alpha=0.5,
                c=cols,
                s=100
            )
        xmin = min(xmin, plottable.points[:,0].min())
        xmax = max(xmax, plottable.points[:,0].max())
        ymin = min(ymin, plottable.points[:,1].min())
        ymax = max(ymax, plottable.points[:,1].max())
        zmin = min(zmin, plottable.points[:,2].min())
        zmax = max(zmax, plottable.points[:,2].max())
        for line in plottable.lines:
            p1 = plottable.points[line[0],:]
            p2 = plottable.points[line[1],:]
            ax.plot([p1[0],p2[0]], [p1[1],p2[1]], [p1[2],p2[2]],c=cols)
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