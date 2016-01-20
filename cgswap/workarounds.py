from cgswap.geometry import PointCloud, TransformationMatrix, plot_3d

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
    def scale_match_vector(self, from_residue, action):
        # We first need the from_vector from the original molecule
        # and the to_vector from the current molecule
        (from_pin1, to_pin1), (from_pin2, to_pin2) = action["pins"].items()
        # From pins are references to the from_residue. These coordinates
        # are correct.
        # To pins are references to the to_residue. These coordinates should be
        # taken from this object instead.
        to_atoms = action["to"]
        from_atoms = action["from"]
        # Convert slice to atom ids
        from_pin1 = from_atoms[from_pin1]
        from_pin2 = from_atoms[from_pin2]
        to_pin1 = to_atoms[to_pin1]
        to_pin2 = to_atoms[to_pin2]

        f1_pos = from_residue.get_mda_atom(from_pin1).position
        f2_pos = from_residue.get_mda_atom(from_pin2).position
        _, t1_pos = self.get_atom_by_id(to_pin1)
        _, t2_pos = self.get_atom_by_id(to_pin2)

        # We want operand_vector to equal operand_vector after scaling
        reference_vector = f1_pos - f2_pos
        operand_vector   = t1_pos - t2_pos
        check = PointCloud(3)
        check.add_points([f1_pos, f2_pos, t1_pos, t2_pos])
        origin = list(t1_pos)
        translation = list(-t1_pos)
        t = TransformationMatrix(3)
        t.translate(*translation)
        s = TransformationMatrix(3)
        s.scale(*list(reference_vector/operand_vector))
        rt = TransformationMatrix(3)
        rt.translate(*origin)
        self.transform(t, subset=to_atoms)
        self.transform(s, subset=to_atoms)
        self.transform(rt, subset=to_atoms)
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
        







