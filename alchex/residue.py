import MDAnalysis as mda
import numpy
from alchex.errors import AtomMismatchException
from alchex.geometry import PointCloud
from copy import deepcopy

class ResidueStructure(object):
    def __init__(self, mdanalysis_residue, residue_parameters):
        self.parameters = residue_parameters
        if type(mdanalysis_residue) is mda.Universe:
            self.mda_object = mdanalysis_residue.select_atoms("all")
        else:
            self.mda_object = mdanalysis_residue
        self.resname = self.mda_object.atoms[0].resname
        self.resid = self.mda_object.atoms[0].resid
        self.verify_residues()
        self.residue_count = 0
    def __getitem__(self, val):
        if type(val) is tuple:
            return [self.mda_object[self.mda_index(x)] for x in val]
        else:
            return self.mda_object[self.mda_index(val[0])]
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
            param_id = self.mda_to_atom_id(atom_id)
            valid = True
            if param_id not in self.parameters.atoms:
                valid = False
            elif self.parameters.atoms[param_id].attrs["atom"] != atom.name:
                valid = False
            if not valid:
                raise AtomMismatchException()
    def mda_to_atom_id(self, mda_index):
        return str(mda_index + 1)
    def mda_index(self, atom_id):
        return int(atom_id) - 1
    def get_mda_atom(self, atom_id):
        return self.mda_object.atoms[self.mda_index(atom_id)]
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

class MultiResidueStructure(ResidueStructure):
    def __init__(self, residue_list):
        start = residue_list[0]
        self.mda_map = {"1."+str(x+1):x for x in range(len(start.mda_object))}
        self.parameters = start.parameters.clone()
        self.parameters.prefix("1.")
        self.mda_object = mda.Merge(start.mda_object).select_atoms("all")
        self.resname = start.resname
        self.resid = -1
        self.resids   = [start.resid]
        self.residue_count = len(residue_list)
        n = 1
        for residue in residue_list[1:]:
            n += 1
            mda_map = {str(n)+"."+str(x+1):x+len(self.mda_object) for x in range(len(residue.mda_object))}
            self.mda_map = dict(mda_map.items() + self.mda_map.items())
            self.resids.append(residue.resid)
            parameters = residue.parameters.clone()
            parameters.prefix(str(n)+".")
            self.parameters.merge(parameters)
            self.mda_object = mda.Merge(self.mda_object, residue.mda_object)
        self.mda_object = self.mda_object.select_atoms("all")
        self.rev_mda_map = {v:k for k, v in self.mda_map.items()}
        self.verify_residues()
    def mda_index(self, atom_id):
        return self.mda_map[atom_id]
    def mda_to_atom_id(self, mda_index):
        return self.rev_mda_map[mda_index]