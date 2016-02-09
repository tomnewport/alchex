# coding: utf-8
from alchex.gromacs_interface import GromacsITPFile, GromacsMDPFile, SimulationContainer, GromacsWrapper
from alchex.exchange_map import ExchangeMap
from alchex.residue import ResidueStructure
from alchex.workarounds import dict_to_gro
from alchex.residue_parameters import ResidueParameters
from alchex.errors import ExchangeMapMissingException
import MDAnalysis as mda
from random import choice
import numpy
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
        self.simulation_containers = {}
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
        #itp        = GromacsITPFile(filename)
        #parameters = itp.read_residue(resname)
        parameters = ResidueParameters(resname)
        parameters.import_itp(filename)
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
    def simulation_container(self, name="alchex_tmp"):
        if name not in self.simulation_containers:
            self.simulation_containers[name] = SimulationContainer(
                name, 
                GromacsWrapper(
                    gromacs_executable=self.gromacs_executable
                    ))
        return self.simulation_containers[name]
    def generate_structure(self, itp_file, resname, repeats=10, additional_itps=[]):
        sims = self.simulation_container("build_"+resname)
        res_params = ResidueParameters(resname)
        res_params.import_itp(itp_file)
        sims.add_file(itp_file)
        topology_string = ""
        for additional_itp in additional_itps:
            topology_string += '#include "{itp_file}"\n'.format(itp_file=path.basename(additional_itp))
            sims.add_file(additional_itp)
        topology_string += '''#include "{itp_file}"
[ system ]
builder
[ molecules ]
{resname}   {repeats}
        '''.format(
            resname=resname, 
            repeats=repeats, 
            itp_file=path.split(itp_file)[1])
        with open(sims.resolve_path("res.top"), "w") as file_handle:
            file_handle.write(topology_string)
        n_atoms = len(res_params.atoms)*repeats
        coordinates = numpy.random.randn(n_atoms,3) * 3 * (n_atoms ** (1./3.))
        coordinates = coordinates - coordinates.min(axis=0)
        box_vector = coordinates.max(axis=0)*1.5
        gro_file = "BUILDER\n{n_atoms}\n".format(n_atoms=n_atoms)
        idx=0
        for rep_id in range(repeats):
            for atom_id in range(len(res_params.atoms)):
                xpos, ypos, zpos = coordinates[idx,:]
                idx += 1
                gro_file += dict_to_gro({
                    "resid":rep_id+1, 
                    "resname":resname,
                    "name":"ANY",
                    "atomid":atom_id+1,
                    "posx":xpos,
                    "posy":ypos,
                    "posz":zpos,
                    "velx":0,
                    "vely":0,
                    "velz":0
                    })+"\n"
        gro_file += "".join([str(x/10)[:9].rjust(10) for x in box_vector])+"\n"
        with open(sims.resolve_path("res.gro"), "w") as file_handle:
            file_handle.write(gro_file)
        em_mdp = '''emstep                   = 0.001
integrator               = steep
emtol                    = 10
nsteps                   = 50000
nstxout                  = 1
        '''
        with open(sims.resolve_path("res.mdp"), "w") as file_handle:
            file_handle.write(em_mdp)
        exitcode, message = sims.gromacs.grompp(kwargs={
            "-f" : "res.mdp",
            "-p" : "res.top",
            "-c" : "res.gro",
            "-o" : "res_em.tpr",
            "-maxwarn":1
            })
        if exitcode != 0:
            print message
            return False
        exitcode, message = sims.gromacs.mdrun(kwargs={
            "-deffnm" : "res_em"
            })
        print message
        if exitcode != 0:
            print message
            return False
        exitcode, message = sims.gromacs.trjconv(kwargs={
            "-f":"res_em.gro",
            "-o":"pbc.gro",
            "-s":"res_em.tpr",
            "-pbc":"mol",
            "_stdin":"0"
        })
        exitcode, message = sims.gromacs.trjconv(kwargs={
            "-f":"res_em.trr",
            "-o":"pbc.xtc",
            "-s":"res_em.tpr",
            "-pbc":"mol",
            "_stdin":"0"
        })

def user_config_path(config):
    return path.join(path.expanduser("~/.alchex/config/"), config,"")

def default_config_path(config):
    return path.join(path.split(__file__)[0], "config", config, "")


def default_configuration():
    folder = path.join(path.split(__file__)[0], "default_configuration")
    #defaultconfig = AlchexConfig(folder=folder, gromacs_executable="/sbcb/packages/opt/Linux_x86_64/gromacs/5.1/bin/gmx_sse")
    defaultconfig = AlchexConfig("cg_default")
    defaultconfig.gromacs_executable = "/sbcb/packages/opt/Linux_x86_64/gromacs/5.1/bin/gmx_avx"
    #defaultconfig.gromacs_executable = "gmx"
    #defaultconfig.generate_structure("CDDG/CDDG2-noconstraints.itp", "CDDG", additional_itps=["data/martini_v2.1-dna.itp"], repeats=1)

    defaultconfig.load_itp_file("data/DLPG.itp", "DLPG")
    defaultconfig.load_itp_file("CDDG/cddg_v3.itp", "CDDG")
    defaultconfig.load_itp_file("data/DVPE.itp", "DVPE")
    defaultconfig.load_itp_file("data/CDL.itp", "CDL0")
    defaultconfig.load_itp_file("data/chol.itp", "CHOL")
    defaultconfig.load_itp_file("data/dppc.itp", "DPPC")
    defaultconfig.load_itp_file("data/pi3p.itp", "PI3P")
    defaultconfig.load_itp_file("data/PODG.itp", "PODG")
    defaultconfig.load_itp_file("data/martini_v2.1.itp", "POPA")
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

    defaultconfig.build_exchange_map(
    from_resname="DPPC",
    from_moltype="DPPC",
    to_resname="PODG",
    to_moltype="PODG",
    exchange_model="martini.lipid")

    defaultconfig.build_exchange_map(
    from_resname="DPPC",
    from_moltype="DPPC",
    to_resname="POPA",
    to_moltype="POPA",
    exchange_model="martini.lipid")

    defaultconfig.build_exchange_map(
    from_resname="DPPC",
    from_moltype="DPPC",
    to_resname="POPG",
    to_moltype="POPG",
    exchange_model="martini.lipid")
    
    defaultconfig.build_exchange_map(
    from_resname="DPPC",
    from_moltype="DPPC",
    to_resname="CDDG",
    to_moltype="CDDG",
    exchange_model="martini.lipid")

    defaultconfig.build_exchange_map(
    from_resname="DPPC",
    from_moltype="DPPC",
    to_resname="DPPC",
    to_moltype="DPPC",
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
    defaultconfig.add_reference_structure("PODG", "data/PODG-em.gro", selection="resname PODG")
    defaultconfig.add_reference_structure("CDDG", "CDDG/cddg_v3.gro")
    defaultconfig.add_reference_structure("POPA", "data/POPA-em.gro")

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

