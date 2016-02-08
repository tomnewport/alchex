# coding: utf-8
from alchex.gromacs_interface import GromacsITPFile, GromacsMDPFile
from alchex.exchange_map import ExchangeMap
from alchex.residue import ResidueStructure
from alchex.errors import ExchangeMapMissingException
import MDAnalysis as mda
from random import choice
from os import path, makedirs
from shutil import rmtree
import json

class AlchexConfig(object):
    def __init__(self, name):
        self.name = name
        # Parameters should be stored as .ITP
        self.parameters = {}
        # Exchange maps can be loaded from .json
        self.exchange_maps = {}
        # Compositions are depracated
        self.compositions = {}
        # Reference structures need to be loaded
        # from .gro files
        self.reference_structures = {}
        # Grompp parameters can be stored 
        # as json
        self.grompp_parameters = {}
    def save(self):
        save_root = user_config_path(self.name)

        if path.exists(save_root):
            rmtree(save_root)

        makedirs(save_root)

        grompp_root = path.join(save_root, "grompp_parameters")
        makedirs(grompp_root)

        for g_name, parameters in self.grompp_parameters.items():
            filename = path.join(grompp_root, g_name+".mdp")
            parameters.to_file(filename)

        parameters_root = path.join(save_root, "topologies")
        makedirs(parameters_root)
        for resname, parameters in self.parameters.items():
            filename = path.join(parameters_root, resname+".itp")
            parameters.export_itp(filename)

        exchange_maps_root = path.join(save_root, "exchange_maps")
        makedirs(exchange_maps_root)
        for from_resname, to_ems in self.exchange_maps.items():
            em_root = path.join(exchange_maps_root, from_resname)
            makedirs(em_root)
            for to_resname, exchange_map in to_ems.items():
                filename = path.join(em_root, from_resname+"-"+to_resname+".json")
                exchange_map.to_file(filename)

        compositions_root = path.join(save_root, "compositions")
        makedirs(compositions_root)

        structures_root = path.join(save_root, "reference_structures")
        makedirs(structures_root)

        for resname, structures in self.reference_structures.items():
            structure_root = path.join(structures_root, resname)
            makedirs(structure_root)
            for idx, structure in enumerate(structures):
                filename = path.join(structure_root, resname+"-"+str(idx)+".gro")
                structure.export_gro(filename)
    def load(self):
        pass
    def load_itp_file(self, filename, resname):
        itp        = GromacsITPFile(filename)
        parameters = itp.read_residue(resname)
        self.parameters[resname] = parameters
    def get_reference_structure(self, resname):
        return choice(self.reference_structures[resname])
    def build_exchange_map(self, from_resname, from_moltype, to_resname, to_moltype, exchange_model, draw=False, **kwargs):
        from_parameters = self.parameters[from_resname]
        to_parameters = self.parameters[to_resname]
        newmap = ExchangeMap()
        newmap.new(
            from_itp=from_parameters, 
            to_itp=to_parameters, 
            to_moltype=to_moltype,
            from_moltype=from_moltype,
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
            raise ExchangeMapMissingException("Cannot find a map to make a {to_resname} from a {from_resname}".format(to_resname=to_resname, from_resname=from_resname))
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
        self.grompp_parameters[name] = params

def user_config_path(config):
    return path.join(path.expanduser("~/.alchex/config/"), config,"")

def default_config_path(config):
    return path.join(path.split(__file__)[0], "config", config, "")


def default_configuration():
    folder = path.join(path.split(__file__)[0], "default_configuration")
    #defaultconfig = AlchexConfig(folder=folder, gromacs_executable="/sbcb/packages/opt/Linux_x86_64/gromacs/5.1/bin/gmx_sse")
    defaultconfig = AlchexConfig("cg_default")
    
    #defaultconfig.gromacs_executable = "/sbcb/packages/opt/Linux_x86_64/gromacs/5.1/bin/gmx_avx"
    defaultconfig.gromacs_executable = "gmx"


    defaultconfig.load_itp_file("data/DLPG.itp", "DLPG")
    defaultconfig.load_itp_file("data/DVPE.itp", "DVPE")
    defaultconfig.load_itp_file("data/CDL.itp", "CDL0")
    defaultconfig.load_itp_file("data/chol.itp", "CHOL")
    defaultconfig.load_itp_file("data/dppc.itp", "DPPC")
    defaultconfig.load_itp_file("data/pi3p.itp", "PI3P")
    defaultconfig.load_itp_file("data/martini_v2.1.itp", "PVPA")
    defaultconfig.load_itp_file("data/martini_v2.1.itp", "POPS")
    defaultconfig.load_itp_file("data/martini_v2.1.itp", "POPE")
    defaultconfig.load_itp_file("data/martini_v2.1.itp", "POPG")
    defaultconfig.load_itp_file("data/popc_patch/popc_ac.itp", "POPC")

    defaultconfig.build_exchange_map(
    from_resname="DPPC", 
    from_moltype="DPPC", 
    to_resname = "CHOL",
    to_moltype = "CHOL", 
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
    from_moltype="POPC",
    to_resname="DPPC",
    to_moltype="DPPC",
    exchange_model="martini.lipid")

    defaultconfig.build_exchange_map(
    from_resname="POPE",
    from_moltype="POPE",
    to_resname="DPPC",
    to_moltype="DPPC",
    exchange_model="martini.lipid")

    defaultconfig.build_exchange_map(
    from_resname="POPE",
    from_moltype="POPE",
    to_resname="POPS",
    to_moltype="POPS",
    exchange_model="martini.lipid")

    defaultconfig.build_exchange_map(
    from_resname="POPE",
    from_moltype="POPE",
    to_resname="DLPG",
    to_moltype="DLPG",
    exchange_model="martini.lipid")

    defaultconfig.build_exchange_map(
    from_resname="POPE",
    from_moltype="POPE",
    to_resname="CDL0",
    to_moltype="CDL0",
    exchange_model="martini.lipid_to_card")

    defaultconfig.build_exchange_map(
    from_resname="DPPC",
    from_moltype="DPPC",
    to_resname="PI3P",
    to_moltype="PI3P",
    exchange_model="martini.lipid")

    defaultconfig.build_exchange_map(
    from_resname="DPPC",
    from_moltype="DPPC",
    to_resname="DPPC",
    to_moltype="DPPC",
    exchange_model="martini.lipid")

    defaultconfig.build_exchange_map(
    from_resname="DPPC",
    from_moltype="DPPC",
    to_resname="CDL0",
    to_moltype="CDL0",
    exchange_model="martini.lipid_to_card")

    defaultconfig.build_exchange_map(
    from_resname="POPC",
    from_moltype="POPC",
    to_resname="CDL0",
    to_moltype="CDL0",
    exchange_model="martini.lipid_to_card")

    defaultconfig.build_exchange_map(
    from_resname="DPPC",
    from_moltype="DPPC",
    to_resname="DLPG",
    to_moltype="DLPG",
    exchange_model="martini.lipid")

    defaultconfig.build_exchange_map(
    from_resname="DPPC",
    from_moltype="DPPC",
    to_resname="POPS",
    to_moltype="POPS",
    exchange_model="martini.lipid")

    defaultconfig.build_exchange_map(
    from_resname="DPPC",
    from_moltype="DPPC",
    to_resname="POPC",
    to_moltype="POPC",
    exchange_model="martini.lipid")

    defaultconfig.build_exchange_map(
    from_resname="DPPC",
    from_moltype="DPPC",
    to_resname="PVPA",
    to_moltype="PVPA",
    exchange_model="martini.lipid")

    defaultconfig.build_exchange_map(
    from_resname="DPPC",
    from_moltype="DPPC",
    to_resname="POPE",
    to_moltype="POPE",
    exchange_model="martini.lipid")

    defaultconfig.build_exchange_map(
    from_resname="POPG",
    from_moltype="POPG",
    to_resname="POPS",
    to_moltype="POPS",
    exchange_model="martini.lipid")

    defaultconfig.build_exchange_map(
    from_resname="POPG",
    from_moltype="POPG",
    to_resname="CDL0",
    to_moltype="CDL0",
    exchange_model="martini.lipid_to_card")

    defaultconfig.build_exchange_map(
    from_resname="POPG",
    from_moltype="POPG",
    to_resname="DPPC",
    to_moltype="DPPC",
    exchange_model="martini.lipid")

    defaultconfig.build_exchange_map(
    from_resname="POPG",
    from_moltype="POPG",
    to_resname="POPG",
    to_moltype="POPG",
    exchange_model="martini.lipid")

    defaultconfig.build_exchange_map(
    from_resname="POPG",
    from_moltype="POPG",
    to_resname="DLPG",
    to_moltype="DLPG",
    exchange_model="martini.lipid")

    defaultconfig.add_reference_structure("DLPG","data/DLPG-em.gro")
    defaultconfig.add_reference_structure("DPPC","data/dppc.pdb", selection="resname DPP")
    defaultconfig.add_reference_structure("CDL0","data/CDL0.gro", selection="resname CDL")
    defaultconfig.add_reference_structure("PI3P","data/pi3p.pdb", selection="resname PI3")
    defaultconfig.add_reference_structure("CHOL","data/chol.pdb", selection="resname CHO")
    defaultconfig.add_reference_structure("PVPA","data/pvpa.pdb", selection="resname PVP")
    defaultconfig.add_reference_structure("POPS","data/pops.pdb", selection="resname POP")
    defaultconfig.add_reference_structure("POPE","data/pope.pdb", selection="resname POP")
    defaultconfig.add_reference_structure("POPC","data/popc.pdb", selection="resname POP")
    defaultconfig.add_reference_structure("POPG", "data/popg.pdb", selection="resname POP")

    defaultconfig.add_grompp_parameters("em", "gromacs_scratch/em.mdp")
    defaultconfig.add_grompp_parameters("alchembed", "alchembed-cg.mdp")

    defaultconfig.save()

    return defaultconfig

f = default_configuration()
'''
f.parameters["POPC"].export_itp("file.itp")
from residue_parameters import ResidueParameters
a = ResidueParameters("POPC")
a.import_itp("martini_v2.1.itp")
a.export_itp("file_out.itp")
'''
f.reference_structures["POPC"][0].export_gro("test_popc.gro")

