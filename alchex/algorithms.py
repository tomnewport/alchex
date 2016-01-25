import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import KDTree
from itertools import combinations
import networkx as nx
from random import random, shuffle


def bruteforce_pair(points, maxiterations=100, d_max = 15, n_tol=0, d_tol=15):

    k = KDTree(points)

    solutions = []

    best_so_far = None

    for a in range(maxiterations):
        g = nx.Graph()
        g.add_nodes_from(range(points.shape[0]))
        g.add_edges_from(k.query_pairs(d_max))
        pairs = []
        i = 0
        min_deg = [1,1]
        while True:
            deg = g.degree()
            di = deg.items()
            shuffle(di)
            min_deg = min(di, key=lambda x : x[1] if x[1] != 0 else 10000)
            i += 1
            n = min_deg[0]
            if min_deg[1] != 0:
                newpairs = g.edges(min_deg[0])
                min_deg_other = 10000
                min_other = -100
                shuffle(newpairs)
                for pair in newpairs:
                    other = pair[1] if pair[1] != min_deg[0] else pair[0]
                    d = deg[other]
                    if d < min_deg_other:
                        min_deg_other = int(d)
                        min_other = other
                newpair = [min_deg[0], min_other]
                pairs.append(newpair)
                g.remove_nodes_from(newpair)
            if i > 2000 or len(g) == 0:
                break
        distance = 0
        if len(pairs) != 0:
            distances = [np.linalg.norm(points[p1] - points[p2]) for p1,p2 in pairs]
            distance = np.max(distances)
        if best_so_far is None or (best_so_far[0] > len(g)) and (best_so_far[1] > distance):
            best_so_far = (len(g), distance, g.nodes(), pairs)
            if len(g) <= n_tol and distance <= d_tol:
                break
    return best_so_far[3], best_so_far[2], distance