from cgswap.geometry import PointCloud, plot_3d

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
        c = to_residue.point_cloud()
        c.transform(tm)
        c.points = c.points[0:1,:]
        plot_3d(from_residue.point_cloud(), to_residue.point_cloud(), c, self.coordinates)
        node_points = PointCloud(3)
        tp = []
        for node in nodes:
            t_index = to_residue.mda_index(node)
            tp.append(to_residue.mda_object.atoms[t_index].position)
            self.ids.append(to_atom_id)
            self.atoms.append(mda_atom_to_dict(to_residue.mda_object.atoms[t_index]))
        node_points.add_points(tp)
        node_points.transform(tm)
        self.coordinates.add_points(node_points.points[:,:3])
    def transform(self, transformation_matrix):
        self.coordinates.transform(transformation_matrix)
