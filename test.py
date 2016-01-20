import networkx as nx
import matplotlib.pyplot as plt

def hamilton(G):
    F = [(G,[G.nodes()[0]])]
    n = G.number_of_nodes()
    while F:
        graph,path = F.pop()
        confs = []
        for node in graph.neighbors(path[-1]):
            conf_p = path[:]
            conf_p.append(node)
            conf_g = nx.Graph(graph)
            conf_g.remove_node(path[-1])
            confs.append((conf_g,conf_p))
        for g,p in confs:
            print(len(p), n)
            if len(p)==n:
                return p
            else:
                F.append((g,p))
    return None

a = nx.Graph()

a.add_nodes_from([1,2,3,4,5,6,7])
a.add_edges_from([[1,2],[2,3],[4,7],[3,4],[4,5],[5,6],[6,7],[7,1]])

g = nx.draw(a)

hamilton(a)

plt.show()


