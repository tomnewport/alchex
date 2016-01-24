from alchex.gromacs import GromacsITPFile
from alchex.exchange_map import ExchangeMap
from alchex.residue import ResidueStructure
import MDAnalysis as mda
from random import choice

class AlchexConfig(object):
    def __init__(self, folder="alchex_configuration"):
        self.folder = folder
        self.parameters = {}
        self.exchange_maps = {}
        self.compositions = {}
        self.reference_structures = {}
    def load_itp_file(self, filename, resname):
        itp        = GromacsITPFile(filename)
        parameters = itp.read_residue(resname)
        self.parameters[resname] = parameters
    def get_reference_structure(self, resname):
        return choice(self.reference_structures[resname])
    def build_exchange_map(self, from_resname, to_resname, exchange_model, draw=False, **kwargs):
        from_parameters = self.parameters[from_resname]
        to_parameters = self.parameters[to_resname]
        newmap = ExchangeMap()
        newmap.new(
            from_itp=from_parameters, 
            to_itp=to_parameters, 
            method=exchange_model, 
            draw=False, 
            **kwargs)
        if from_resname not in self.exchange_maps:
            self.exchange_maps[from_resname] = {}
        self.exchange_maps[from_resname][to_resname] = newmap
        print(newmap)
    def add_composition(self, name, **resname_fractions):
        n_val = sum(resname_fractions.values()) * 1.0
        resname_fractions = {k: v/n_val for k, v in resname_fractions.items()}
        self.compositions[name] = resname_fractions
    def add_reference_structure(self, resname, structure_file):
        ref_universe  = mda.Universe(structure_file)
        if resname not in self.reference_structures:
            self.reference_structures[resname] = []
        for residue in ref_universe.select_atoms("resname "+resname).residues:
            structure = ResidueStructure(mda.Merge(residue), self.parameters[resname])
            self.reference_structures[resname].append(structure)










