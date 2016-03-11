from knowledge import martini_atom_similarity
from itertools import product
import matplotlib.pyplot as plt
import networkx as nx
from copy import deepcopy
from alchex.knowledge import ITP_FIELDS
import re

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
    def to_itp_line(self):
        fields = []
        for field_name in ITP_FIELDS["atoms"]:
            if field_name in self.attrs:
                fields.append(self.attrs[field_name])
            else:
                fields.append("")
        return "\t".join(fields)



class RPBond(object):
    table_name = "bonds"
    def __init__(self, i_atom, j_atom, attrs):
        i_atom.i_bonded.append(self)
        j_atom.j_bonded.append(self)
        self.attrs = attrs
        self.i_atom = i_atom
        self.j_atom = j_atom
    def to_itp_line(self):
        fields = []
        for field_name in ITP_FIELDS[self.table_name]:
            if field_name in self.attrs:
                fields.append(self.attrs[field_name])
            else:
                fields.append("")
        return "\t".join(fields)

class RPConstraint(RPBond):
    table_name = "constraints"

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
    def to_itp_line(self):
        fields = []
        for field_name in ITP_FIELDS["angles"]:
            if field_name in self.attrs:
                fields.append(str(self.attrs[field_name]))
            else:
                fields.append("")
        return "\t".join(fields)

class RPDihedral(object):
    def __init__(self, attrs):
        self.attrs = attrs

class ResidueParameters(object):
    def __init__(self, resname):
        self.resname = resname
        self.nrexcl = "1"
        self.atoms = {}
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.constraints = []
        self._graph = None
        self._atom_to_id_dict = None
    def export_itp(self, filename):
        itp_text = "[ moleculetype ]\n"
        itp_text += self.resname+" "+self.nrexcl+"\n\n"
        itp_text += "[ atoms ]\n"
        for atom_id, atom in sorted(self.atoms.items(), key=lambda x : int(x[0])):
            itp_text +=  atom.to_itp_line()+"\n"
        itp_text += "\n[ bonds ]\n"
        for bond in self.bonds:
            itp_text +=  bond.to_itp_line()+"\n"
        itp_text += "\n[ angles ]\n"
        for angle in self.angles:
            itp_text +=  angle.to_itp_line()+"\n"
        with open(filename, "w") as file_handle:
            file_handle.write(itp_text)
    def import_itp(self, itp_file):
        with open(itp_file, "r") as file_handle:
            table_head = re.compile(r'^\s*\[\s*(\w+)\s*\]\s*$')
            section_started = False
            current_table_name = ""
            for line in file_handle:
                table_head_search = table_head.search(line)
                if table_head_search is not None:
                    current_table_name = table_head_search.groups(1)[0]
                else:
                    if ";" in line:
                        line = line.split(";")[0]
                    line = line.strip()
                    if line != "" and line[0] != "#":
                        parts = line.split()
                        if current_table_name == "moleculetype":
                            if section_started:
                                return
                            if not section_started:
                                if parts[0] == self.resname:
                                    section_started = True
                                    self.nrexcl = parts[1]
                        elif section_started:
                            if current_table_name == "atoms":
                                atom_attrs = {k:v for k, v in zip(ITP_FIELDS["atoms"], parts)}
                                self.add_atom(atom_attrs["id"],atom_attrs)
                            elif current_table_name == "bonds":
                                bond_attrs = {k:v for k, v in zip(ITP_FIELDS["bonds"], parts)}
                                self.add_bond(bond_attrs["i"], bond_attrs["j"],bond_attrs)
                            elif current_table_name == "angles":
                                angle_attrs = {k:v for k, v in zip(ITP_FIELDS["angles"], parts)}
                                self.add_angle(angle_attrs["i"], angle_attrs["j"],angle_attrs["k"],angle_attrs)
                            elif current_table_name == "constraints":
                                constraint_attrs = {k:v for k, v in zip(ITP_FIELDS["constraints"], parts)}
                                self.add_constraint(constraint_attrs["i"], constraint_attrs["j"],constraint_attrs)
                            elif current_table_name == "dihedrals":
                                dh_attrs = {k:v for k, v in zip(ITP_FIELDS["dihedrals"], parts)}
                                self.add_dihedral(dh_attrs["i"], dh_attrs["j"],dh_attrs["k"],dh_attrs["l"], dh_attrs)

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
    def add_constraint(self, i, j, attributes):
        i_atom = self.atoms[i]
        j_atom = self.atoms[j]
        self.constraints.append(RPConstraint(i_atom, j_atom, attributes))
    def add_angle(self, i, j, k, attributes):
        i_atom = self.atoms[i]
        j_atom = self.atoms[j]
        k_atom = self.atoms[k]
        self.angles.append(RPAngle(i_atom, j_atom, k_atom, attributes))
    def add_dihedral(self, i,j,k,l,attributes):
        self.dihedrals.append(RPDihedral(attributes))
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