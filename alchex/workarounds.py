from alchex.geometry import PointCloud, TransformationMatrix, plot_3d
from scipy.spatial.distance import cdist
import numpy
import networkx as nx
from itertools import combinations
from copy import deepcopy

def mda_atom_to_dict(mda_atom):
    return {
        "name" : mda_atom.name
    }

class WAEditableResidue(object):
    def __init__(self, resname, resid):
        self.resname = resname
        self.resid = resid
        self.ids = []
        self.atoms = []
        self.coordinates = PointCloud(3)
    def import_mdanalysis_atoms(self, residue, id_mapping="direct"):
        if id_mapping == "direct":
            for atom_id, atom in enumerate(residue.mda_object.atoms):
                self.ids.append(str(atom_id + 1))
                self.atoms.append(mda_atom_to_dict(atom))
            self.coordinates.add_points(residue.mda_object.coordinates())
            #plot_3d(self.coordinates)
    def overlay(self, from_residue, to_residue, mapping):
        points = []
        for from_atom_id, to_atom_id in mapping.items():
            f_index = from_residue.mda_index(from_atom_id)
            t_index = to_residue.mda_index(to_atom_id)
            # We want to add an atom:
            #    with the index and name of to_residue
            #    with the position of from_residue
            self.ids.append(to_atom_id)
            self.atoms.append(mda_atom_to_dict(to_residue.mda_object.atoms[t_index]))
            points.append(from_residue.mda_object.atoms[f_index].position)
        self.coordinates.add_points(points)
    def align_fragment(self, from_residue, to_residue, params):
        nodes = params["nodes"]
        centroid_weighting = params["centroid_weighting"]
        reference = params["reference"]
        # Nodes refers to atoms of the to_residue
        # Reference contains a mapping of from, to
        # We need to:
        #    1. compute rotation to map reference pointclouds 
        #       onto each other - SVD alignment
        #    2. compute translation to map weighted centroids
        #       onto each other - weighted means
        from_points = PointCloud(3)
        to_points = PointCloud(3)
        fp = []
        tp = []
        for from_atom_id, to_atom_id in reference:
            f_index = from_residue.mda_index(from_atom_id)
            t_index = to_residue.mda_index(to_atom_id)
            tp.append(to_residue.mda_object.atoms[t_index].position)
            fp.append(from_residue.mda_object.atoms[f_index].position)
        from_points.add_points(fp)
        to_points.add_points(tp)
        tm, rmse, aligned = from_points.paired_3d_align(to_points, centroid_weighting=centroid_weighting, inv=False)
        # tm now corrects rotation with centroid weighting function
        #c = to_residue.point_cloud()
        #c.transform(tm)
        #c.points = c.points[0:1,:]
        #plot_3d(from_residue.point_cloud(), to_residue.point_cloud(), c, self.coordinates)
        node_points = PointCloud(3)
        tp = []
        for node in nodes:
            t_index = to_residue.mda_index(node)
            tp.append(to_residue.mda_object.atoms[t_index].position)
            self.ids.append(node)
            self.atoms.append(mda_atom_to_dict(to_residue.mda_object.atoms[t_index]))
        node_points.add_points(tp)
        node_points.transform(tm)
        self.coordinates.add_points(node_points.points[:,:3])
    def get_atom_by_id(self, atom_id):
        index = self.ids.index(atom_id)
        return self.atoms[index], self.coordinates.points[index,:3]
    def transform(self, transformation_matrix, subset = None):
        if subset is None:
            self.coordinates.transform(transformation_matrix)
        else:
            tcoords = self.coordinates.clone()
            tcoords.transform(transformation_matrix)
            for atom_id in subset:
                idx = self.ids.index(atom_id)
                self.coordinates.points[idx, :] = tcoords.points[idx, :]
    def stretch_interpolate(self, from_residue, action):
        # First extract the from_residue points:
        from_points = []
        for fpid in action["from"]:
            from_points.append(from_residue.get_mda_atom(fpid).position)
        from_pc = PointCloud(3)
        from_pc.add_points(from_points)
        new_pc = from_pc.interpolate_1d_list(numpy.linspace(0,1,len(action["to"])))
        for idx, new_atom_id in enumerate(action["to"]):
            self.coordinates.points[self.ids.index(new_atom_id),:] = new_pc.points[idx,:]
    def build_bridge(self, from_residue, to_residue, action):
        from_atom = action["from"]
        to_atom = action["to"]
        middle_atoms = action["atoms"][1:-1]
        _, from_atom_position = self.get_atom_by_id(from_atom)
        _, to_atom_position = self.get_atom_by_id(to_atom)
        ends = PointCloud(3)
        ends.add_points([from_atom_position, to_atom_position])
        points = []
        for seq_id, atom in enumerate(middle_atoms):
            coordinate_1d = (seq_id + 1.0) / (len(action["atoms"]) - 1.0)
            self.ids.append(atom)
            self.atoms.append(mda_atom_to_dict(to_residue.get_mda_atom(atom)))
            point = ends.centroid(centroid_weighting=[coordinate_1d,1-coordinate_1d])[:3]
            points.append(point)
        self.coordinates.add_points(points)
    def as_pdb(self):
        lines = []
        for atom_id in self.ids:
            print atom_id
            atom, coordinates = self.get_atom_by_id(atom_id)
            line = "ATOM    {id} {name} {resname} {resid}      {xpos}  {ypos}  {zpos}  1.00  0.00           N"
            lines.append(line.format(
                id=str(atom_id).rjust(3), 
                name=atom["name"].rjust(4), 
                resname=self.resname.rjust(4),
                resid = str(self.resid).rjust(4), 
                xpos=str(coordinates[0])[:6].ljust(6),
                ypos=str(coordinates[1])[:6].ljust(6),
                zpos=str(coordinates[2])[:6].ljust(6)))
        return lines
    def as_gro(self):
        lines = {}
        for atom_id in self.ids:
            atom, coordinates = self.get_atom_by_id(atom_id)
            line = ""
            #residue number (5 positions, integer)
            line += str(self.resid).rjust(5)
            #residue name (5 characters)
            line += str(self.resname).ljust(5)
            #atom name (5 characters)
            line += str(atom["name"]).rjust(5)
            #atom number (5 positions, integer)
            line += str(atom_id).rjust(5)
            #position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places)
            line += "".join(["{:8.4f}".format(x) for x in coordinates/10])
            #velocity (in nm/ps (or km/s), x y z in 3 columns, each 8 positions with 4 decimal places)
            line += "".join(["{:8.4f}".format(x) for x in [0,0,0]])
            lines[int(atom_id)] = line
        lines = sorted(lines.items(), key=lambda x : x[0])
        return [x[1] for x in lines]

def gro_to_dict(groline):
    #residue number (5 positions, integer)
    resid = groline[:5].strip()
    #residue name (5 characters)
    resname = groline[5:10].strip()
    #atom name (5 characters)
    atomname = groline[10:15].strip()
    #atom number (5 positions, integer)
    atomid = groline[15:20].strip()
    #position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places)
    posx = float(groline[20:28])
    posy = float(groline[28:36])
    posz = float(groline[36:44])
    #velocity (in nm/ps (or km/s), x y z in 3 columns, each 8 positions with 4 decimal places)
    velx = float(groline[44:52])
    vely = float(groline[52:60])
    velz = float(groline[60:68])
    return {
            "resid"  :resid, 
            "resname":resname, 
            "name"   :atomname, 
            "atomid" :atomid, 
            "posx"   :posx, 
            "posy"   :posy, 
            "posz"   :posz, 
            "velx"   :velx, 
            "vely"   :vely, 
            "velz"   :velz}

def dict_to_gro(dictionary):
    line = ""
    #residue number (5 positions, integer)
    line += str(dictionary["resid"]).rjust(5)
    #residue name (5 characters)
    line += str(dictionary["resname"]).ljust(5)
    #atom name (5 characters)
    line += str(dictionary["name"]).rjust(5)
    #atom number (5 positions, integer)
    line += str(dictionary["atomid"]).rjust(5)
    #position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places)
    line += "".join(["{:8.4f}".format(x) for x in [dictionary["posx"], dictionary["posy"], dictionary["posz"]]])
    #velocity (in nm/ps (or km/s), x y z in 3 columns, each 8 positions with 4 decimal places)
    line += "".join(["{:8.4f}".format(x) for x in [dictionary["velx"], dictionary["vely"], dictionary["velz"]]])
    return line

class WAEditableGrofile(object):
    def __init__(self):
        self.sysname = "Edited Gro File"
        self.residues = []
        self.box_vector = [0,0,0]
    def top_molecules(self):
        molecules = []
        for residue in self.residues:
            if len(molecules) == 0 or molecules[-1][0] != residue.resname:
                molecules.append((residue.resname, 1))
            else:
                molecules[-1] = (residue.resname, molecules[-1][1] + 1)
        return molecules
    def clashgraph(self, distance_tolerance=1):
        graph = nx.Graph()
        graph.add_nodes_from(range(len(self.residues)))
        for n1, n2 in combinations(graph, 2):
            p1 = self.residues[n1].coordinates
            p2 = self.residues[n2].coordinates
            distance = cdist(p1.points, p2.points).min()
            if distance < distance_tolerance:
                graph.add_edge(n1, n2)
        return graph
    def declash(self, distance_tolerance=1):
        graphs = []
        clashgraph = self.clashgraph(distance_tolerance)
        for node in clashgraph:
            found = False
            edges = set([x[1] for x in clashgraph.edges(node)])
            for graph, clashes in graphs:
                if node not in clashes:
                    found = True
                    graph.append(node)
                    clashes.update(edges)
                    break
            if not found:
                graphs.append(([node],edges))
        declashed = []
        for idx, (nodes, _) in enumerate(graphs):
            new_grofile = WAEditableGrofile()
            new_grofile.sysname = self.sysname + " (declash "+str(idx)+")"
            new_grofile.box_vector = self.box_vector
            for node in nodes:
                new_grofile.residues.append(self.residues[node])
            declashed.append(new_grofile)
        return declashed
    def from_file(self, filename):
        with open(filename, "r") as fh:
            lines = [x for x in fh.read().split("\n") if x.strip() != ""]
        self.sysname = lines[0].strip()
        for line in lines[2:-1]:
            atom_data = gro_to_dict(line)
            coordinates = [[
                10*float(atom_data["posx"]), 
                10*float(atom_data["posy"]), 
                10*float(atom_data["posz"])
                ]]
            if atom_data["resid"] not in self.residues:
                self.residues.append(WAEditableResidue(resid=atom_data["resid"], resname=atom_data["resname"]))
                res_start = int(atom_data["atomid"])
            self.residues[-1].coordinates.add_points(coordinates)
            atom_id = str(1 + int(atom_data["atomid"]) - res_start)
            self.residues[-1].ids.append(atom_id)
            self.residues[-1].atoms.append(atom_data)
        self.box_vector = [10*float(x) for x in lines[-1].split()]
    def list_residues(self):
        rlist = []
        for residue in self.residues:
            if len(rlist) == 0 or residue.resname != rlist[-1][0]:
                rlist.append([residue.resname, 0])
            rlist[-1][1] += 1
        return rlist
    def to_file(self, filename):
        atom_count = sum([len(x.ids) for x in self.residues])
        with open(filename, "w") as file_handle:
            file_handle.write(self.sysname + "\n")
            file_handle.write(str(atom_count).rjust(5)+ "\n")
            atom_id = 0
            for residue in self.residues:
                for idx, atom_data in enumerate(residue.atoms):
                    atom_id += 1
                    atom = {}
                    position = residue.coordinates.points[idx,:3]
                    atom["posx"], atom["posy"], atom["posz"] = position/10
                    atom["atomid"]  = str(atom_id)
                    atom["resid"]   = residue.resid
                    atom["resname"] = residue.resname
                    atom["name"]    = atom_data["name"]
                    atom["velx"]    = atom_data.get("velx", 0)
                    atom["vely"]    = atom_data.get("vely", 0)
                    atom["velz"]    = atom_data.get("velz", 0)
                    file_handle.write(dict_to_gro(atom) + "\n")
            file_handle.write("".join([str(x/10).rjust(10) for x in self.box_vector])+"\n")