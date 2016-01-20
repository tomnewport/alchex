import MDAnalysis as mda
import numpy
from cgswap.errors import AtomMismatchException
from cgswap.geometry import PointCloud

class ResidueStructure(object):
    def __init__(self, mdanalysis_residue, residue_parameters):
        self.parameters = residue_parameters
        self.mda_object = mdanalysis_residue
        self.resname = self.mda_object.atoms[0].resname
        self.resid = self.mda_object.atoms[0].resid
        self.verify_residues()
    def __getitem__(self, val):
        if type(val) is tuple:
            return [self.mda_object[int(x)-1] for x in val]
        else:
            return self.mda_object[int(val[0])-1]
    def transform(self,matrix):
        p = self.point_cloud()
        p.transform(matrix)
        for atom_id, atom in enumerate(self.mda_object.atoms):
            atom.position = p.points[atom_id, :3]
    def point_cloud(self):
        p = PointCloud(3)
        p.add_points(self.mda_object.coordinates())
        return p
    def position(self, atom_id):
        if type(atom_id) is str:
            return self[atom_id].position
        else:
            return numpy.array([self[x].position for x in atom_id])
    def verify_residues(self):
        for atom_id, atom in enumerate(self.mda_object.atoms):
            param_id = str(atom_id + 1)
            valid = True
            if param_id not in self.parameters.atoms:
                valid = False
            elif self.parameters.atoms[param_id].attrs["atom"] != atom.name:
                valid = False
            if not valid:
                print(atom.name)
                raise AtomMismatchException()
    def mda_index(self, atom_id):
        return int(atom_id) - 1
    def change_residue_name(self, resname):
        self.resname = resname
        for atom in self.mda_object.atoms:
            atom.resname = resname
    def change_residue_id(self, resid):
        self.resid = resid
        for atom in self.mda_object.atoms:
            atom.resid = resid
    def clone(self):
        return ResidueStructure(mda.Merge(self.mda_object).select_atoms("all"), self.parameters)