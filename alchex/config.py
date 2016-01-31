from alchex.gromacs_interface import GromacsITPFile, GromacsMDPFile
from alchex.exchange_map import ExchangeMap
from alchex.residue import ResidueStructure
from alchex.errors import ExchangeMapMissingException
import MDAnalysis as mda
from random import choice
from os import path

class AlchexConfig(object):
    def __init__(self, 
        folder="alchex_configuration", 
        gromacs_executable="gmx"):
        self.folder = folder
        if path.exists(path.join(folder, "alchex.json")):
            raise NotImplementedError()
        else:
            self.gromacs_executable = gromacs_executable
            self.parameters = {}
            self.exchange_maps = {}
            self.compositions = {}
            self.reference_structures = {}
            self.grompp_parameters = {}
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
    def get_exchange_map(self, from_resname, to_resname):
        if from_resname in self.exchange_maps and to_resname in self.exchange_maps[from_resname]:
            return self.exchange_maps[from_resname][to_resname]
        else:
            raise ExchangeMapMissingException()
    def add_composition(self, name, **resname_fractions):
        n_val = sum(resname_fractions.values()) * 1.0
        resname_fractions = {k: v/n_val for k, v in resname_fractions.items()}
        self.compositions[name] = resname_fractions
    def add_reference_structure(self, resname, structure_file, selection=None):
        ref_universe  = mda.Universe(structure_file)
        if resname not in self.reference_structures:
            self.reference_structures[resname] = []
        if selection is None:
            selection = "resname " + resname
        for residue in ref_universe.select_atoms(selection).residues:
            residue.name = resname
            for atom in residue.atoms:
                atom.resname = resname
            structure = ResidueStructure(mda.Merge(residue), self.parameters[resname])
            self.reference_structures[resname].append(structure)
    def add_grompp_parameters(self, name, mdp_file):
        params = GromacsMDPFile()
        params.from_file(mdp_file)
        self.grompp_parameters[name] = params.attrs

def default_configuration():
    folder = path.join(path.split(__file__)[0], "default_configuration")
    #defaultconfig = AlchexConfig(folder=folder, gromacs_executable="/sbcb/packages/opt/Linux_x86_64/gromacs/5.1/bin/gmx_sse")
    defaultconfig = AlchexConfig(folder=folder, gromacs_executable="gmx")



    defaultconfig.load_itp_file("data/DLPG.itp", "DLPG")
    defaultconfig.load_itp_file("data/DVPE.itp", "DVPE")
    defaultconfig.load_itp_file("data/CDL.itp", "CDL0")
    defaultconfig.load_itp_file("data/chol.itp", "CHOL")
    defaultconfig.load_itp_file("data/dppc.itp", "DPPC")
    defaultconfig.load_itp_file("data/dppc.itp", "DPPC")
    defaultconfig.load_itp_file("data/pi3p.itp", "PI3P")
    defaultconfig.load_itp_file("data/popc_patch/popc_ac.itp", "POPC")

    defaultconfig.build_exchange_map(
    from_resname="DPPC", 
    to_resname = "CHOL", 
    exchange_model="martini.static_planar_alignment", 
    draw=False, 
    clusters=[
        [
            ["PO4", "NC3"],
            ["ROH"],
            1
        ],
        [
            ["C1B", "C2A"],
            ["R2", "R4"],
            1
        ],
        [
            ["C4B", "C4A"],
            ["C2"],
            1
        ]
    ])

    defaultconfig.build_exchange_map(
    from_resname="POPC",
    to_resname="DPPC",
    exchange_model="martini.lipid")

    defaultconfig.build_exchange_map(
    from_resname="DPPC",
    to_resname="PI3P",
    exchange_model="martini.lipid")

    defaultconfig.build_exchange_map(
    from_resname="POPC",
    to_resname="CDL0",
    exchange_model="martini.lipid_to_card")

    defaultconfig.build_exchange_map(
    from_resname="POPC",
    to_resname="DLPG",
    exchange_model="martini.lipid")

    #defaultconfig.add_reference_structure("DLPG","data/DLPG-em.gro")
    defaultconfig.add_reference_structure("DPPC","data/dppc.pdb", selection="resname DPP")
    defaultconfig.add_reference_structure("CDL0","data/CDL0.gro", selection="resname CDL")
    defaultconfig.add_reference_structure("PI3P","data/pi3p.pdb", selection="resname PI3")

    defaultconfig.add_grompp_parameters("em", "gromacs_scratch/em.mdp")
    defaultconfig.add_grompp_parameters("alchembed", "alchembed-cg.mdp")

    return defaultconfig