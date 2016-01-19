import MDAnalysis as mda
from cgswap.errors import AtomMismatchException

class ResidueStructure(object):
    def __init__(self, mdanalysis_residue, residue_parameters):
        self.parameters = residue_parameters
        self.names = mdanalysis_residue.names
        self.coordinates = mdanalysis_residue.coordinates()
        self.mda_object = mdanalysis_residue
        self.resname = self.mda_object.atoms[0].resname
        self.resid = self.mda_object.atoms[0].resid
        self.number_residues()
    def number_residues(self):
        self.residue_numbering = {}
        for name in self.names:
            try:
                ids = self.parameters.atoms_to_ids(name)
                if len(ids) > 1:
                    raise NotImplementedError
                else:
                    self.residue_numbering[name] = ids[0]
            except:
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