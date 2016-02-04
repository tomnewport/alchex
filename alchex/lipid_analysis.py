#!/usr/bin/env python
# -*- coding: utf-8 -*-

import networkx as nx
from alchex.geometry import PointCloud
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
import numpy

def find_bilayer_leaflets(universe, headgroups="name P*", axis=[0,0,1]):
	points = universe.select_atoms(headgroups)
	coordinates = points.coordinates()
	graph = nx.Graph()
	graph.add_nodes_from(range(len(points)))
	distances = squareform(pdist(coordinates))
	numpy.fill_diagonal(distances, 100000000)
	min_dist = distances.min(axis=1).max()
	for n1, n2_dists in enumerate(distances):
		ranked = numpy.argsort(n2_dists)
		closest = ranked[0]
		for n2 in numpy.argwhere(n2_dists <= min_dist):
			graph.add_edge(n1, n2[0])
	subgraphs = list(nx.connected_component_subgraphs(graph))
	subgraphs = sorted(subgraphs, key=len)
	subgraphs = subgraphs[-2:]
	subgraphs=sorted(subgraphs, key= lambda x : numpy.dot(axis, coordinates.take(x, axis=0).mean(axis=0)))
	# lower, upper
	leaflets = [
		[points[x].resid for x in sg.nodes()] for sg in subgraphs
	]
	return leaflets
		