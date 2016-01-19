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
    def overlay(self, to_residue, from_atom, to_atom):
        # When called on new_residue, this method
        # determines the name of the to_atom in to_residue
        # and renames from_atom of self.mda_object accordingly
        from_atom_id = int(from_atom) - 1
        to_atom = to_residue.parameters.atoms[to_atom]
        self.mda_object[from_atom_id].name = to_atom.attrs["atom"]
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
        self.residue_numbering = {}
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