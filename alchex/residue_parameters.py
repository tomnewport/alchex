from knowledge import martini_atom_similarity
from itertools import product
import matplotlib.pyplot as plt
import networkx as nx
from copy import deepcopy

class RPAtom(object):
    def __init__(self, attrs):
        self.attrs = attrs
        self.i_bonded = []
        self.j_bonded = []
        self.i_angled = []
        self.j_angled = []
        self.k_angled = []
    def prefix(self, prefix):
        self.attrs["id"] = prefix + self.attrs["id"]
    def bonded(self):
        ba = []
        for bond in self.i_bonded:
            ba.append(bond.j_atom)
        for bond in self.j_bonded:
            ba.append(bond.i_atom)
        return ba
    def bonded_ids(self):
        return [x.attrs["id"] for x in self.bonded()]

class RPBond(object):
    def __init__(self, i_atom, j_atom, attrs):
        i_atom.i_bonded.append(self)
        j_atom.j_bonded.append(self)
        self.attrs = attrs
        self.i_atom = i_atom
        self.j_atom = j_atom

class RPAngle(object):
    def __init__(self, i_atom, j_atom, k_atom, attrs):
        i_atom.i_angled.append(self)
        j_atom.j_angled.append(self)
        k_atom.k_angled.append(self)
        self.attrs = attrs
        self.i_atom = i_atom
        self.j_atom = j_atom
        self.k_atom = k_atom
        self._graph = None

class ResidueParameters(object):
    def __init__(self, resname):
        self.resname = resname
        self.atoms = {}
        self.bonds = []
        self.angles = []
        self._graph = None
        self._atom_to_id_dict = None
    def clone(self):
        return deepcopy(self)
    def merge(self, other):
        self.atoms = {k:v for k, v in self.atoms.items() + other.atoms.items()}
        self.bonds += other.bonds
        self.angles += other.angles
        self._graph = None
        self._atom_to_id_dict = None
    def prefix(self, prefix):
        self.atoms = {prefix+k:v for k,v in self.atoms.items()}
        for atom in self.atoms.values():
            atom.prefix(prefix)
    def add_atom(self, id, attributes):
        self.atoms[id] = RPAtom(attributes)
    def add_bond(self, i, j, attributes):
        i_atom = self.atoms[i]
        j_atom = self.atoms[j]
        self.bonds.append(RPBond(i_atom, j_atom, attributes))
    def add_angle(self, i, j, k, attributes):
        i_atom = self.atoms[i]
        j_atom = self.atoms[j]
        k_atom = self.atoms[k]
        self.angles.append(RPAngle(i_atom, j_atom, k_atom, attributes))
    def atoms_to_ids(self, atoms):
        if self._atom_to_id_dict is None:
            self._atom_to_id_dict = {}
            for atom in self.atoms.values():
                if atom.attrs["atom"] not in self._atom_to_id_dict:
                    self._atom_to_id_dict[atom.attrs["atom"]] = []
                self._atom_to_id_dict[atom.attrs["atom"]].append(atom.attrs["id"])
        if type(atoms) is list:
            return [self._atom_to_id_dict[x] for x in atoms]
        else:
            return self._atom_to_id_dict[atoms]
    def ids_to_atoms(self, ids):
        if type(ids) is list:
            return [self.atoms[i].attrs["atom"] for i in ids]
        else:
            return self.atoms[ids].attrs["atom"]
    def bonded_ids(self, atom_id):
        return self.atoms[atom_id].bonded_ids()
    def graph(self):
        if self._graph is None:
            g = nx.Graph()
            atoms = self.atoms.keys()
            g.add_nodes_from(self.atoms.keys())
            for atom in atoms:
                for other in self.bonded_ids(atom):
                    g.add_edge(atom, other)
            self._graph = g
        return self._graph
    def mcs(self, other, threshold=0.9):
        # Start with dictionaries of atoms and what type they're connected to:
        from_atoms = {
            k:(
                atom, 
                atom.bonded()
            ) for k, atom in self.atoms.items()
        }
        to_atoms = {
            k:(
                atom, 
                atom.bonded()
            ) for k, atom in other.atoms.items()
        }
        # These will hold our maps as they're grown
        grow_maps = list()
        done_maps = list()
        # Find common subgraphs to start from (of size >= 3)
        # Iterate over pairs of from atoms and to atoms:
        for from_atom_id, (from_atom, from_connected_atoms) in from_atoms.items():
            for to_atom_id, (to_atom, to_connected_atoms) in to_atoms.items():
                if martini_atom_similarity(from_atom, to_atom) > 0.5:
                    # Atom types are similar
                    candidate_map = {from_atom_id:to_atom_id}
                    connected_similarity = 0
                    # Sum up compatible connected nodes
                    for fc_atom, tc_atom in product(from_connected_atoms, to_connected_atoms):
                        similarity = martini_atom_similarity(fc_atom, tc_atom)
                        connected_similarity += similarity                           
                    if connected_similarity > 1:
                        # If total is greater than 1, add it to the queue
                        grow_maps.append(candidate_map)
        iters = 0
        while len(grow_maps) > 0 and iters < 1000:
            # Prevent possible while loop explosion
            iters += 1
            # Pop next map from front of queue
            a_b_map = grow_maps.pop(0)
            grown = False
            # Compute reverse map
            b_a_map = {v:k for k, v in a_b_map.items()}
            # Try to grow from each atom
            for a_atom, b_atom in a_b_map.items():
                # Find bonded atoms not already in the map
                a_bonded_atoms = [atom for atom in self.atoms[a_atom].bonded() 
                            if atom.attrs["id"] not in a_b_map]
                b_bonded_atoms = [atom for atom in other.atoms[b_atom].bonded() 
                            if atom.attrs["id"] not in b_a_map]
                # Iterate over all pairs of atoms bonded to atoms in our two molecules
                for a_bonded_atom in a_bonded_atoms:
                    for b_bonded_atom in b_bonded_atoms:
                        # Check if similar
                        if martini_atom_similarity(a_bonded_atom, b_bonded_atom) > 0.5:
                            # Copy the dictionary
                            new_map = {x:y for x, y in a_b_map.items()}
                            # Add to the copied dictionary
                            new_map[a_bonded_atom.attrs["id"]] = b_bonded_atom.attrs["id"]
                            # Add unless it's already there
                            if new_map not in grow_maps:
                                grow_maps.append(new_map)
                            grown = True
            # Graph could not be grown further, add to completed list
            if not grown and a_b_map not in done_maps:
                done_maps.append(a_b_map)
                # Control list size - remove smaller graphs
                if len(done_maps) > 1000:
                    done_maps = sorted(done_maps, key=lambda x : -len(x))[:1000]
        sorted_maps = sorted(done_maps, key=lambda x : -len(x))
        longest = sorted_maps[0]

        return [candidate for candidate in sorted_maps if len(candidate) >= len(longest) * threshold]
    def centre_node(self):
        return [x[0] for x in sorted(nx.eccentricity(self.graph()).items(), key=lambda x : x[1])][0][0]
    def path_to_centre(self, node):
        return nx.shortest_path(self.graph(), source=node, target=self.centre_node())
    def eccentricity(self, **kwargs):
        return nx.eccentricity(self.graph(), **kwargs)
    def find_lipid_tails(self):
        tail_components = ["C1","C3"]
        graph = self.graph()
        sgs = nx.connected_component_subgraphs(graph.subgraph(
                [x for x,a in self.atoms.items() if a.attrs["type"] in tail_components]
            ))
        sgs = [x.nodes() for x in sgs]
        subgraphs = []
        for subgraph in sgs:
            subgraphs.append(sorted(subgraph, key=lambda x : len(self.path_to_centre(x))))
        # Returns list of lists of tail nodes, with element 0 closest to molecule's centre
        return subgraphs
    def draw(self, default_colour="white", **colours):
        plt.figure(figsize=(12,12))
        plt.axis('equal')
        g   = self.graph()
        pos = nx.spring_layout(g)
        colour_specified = set([x for y in colours.values() for x in y])
        colours[default_colour] = set(self.atoms.keys()) - colour_specified
        for colour_name, atom_set in colours.items():
            nx.draw_networkx_nodes(g,pos,
                       nodelist=atom_set,
                       node_color=colour_name,
                       node_size=600,
                   alpha=1)
        nx.draw_networkx_edges(g, pos, width=15,alpha=0.5,edge_color='k')
        nx.draw_networkx_labels(g,pos,{k: v.attrs["atom"] for k,v in self.atoms.items()},font_size=8)
        plt.show()