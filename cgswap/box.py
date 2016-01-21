import MDAnalysis as mda
from random import shuffle
from cgswap.residue import ResidueStructure

def gro_to_dict(groline):
    #residue number (5 positions, integer)
    resid = groline[:5]
    #residue name (5 characters)
    resname = groline[5:10]
    #atom name (5 characters)
    atomname = groline[10:15]
    #atom number (5 positions, integer)
    atomid = groline[15:20]
    #position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places)
    posx = groline[20:28]
    posy = groline[28:36]
    posz = groline[36:44]
    #velocity (in nm/ps (or km/s), x y z in 3 columns, each 8 positions with 4 decimal places)
    velx = groline[44:52]
    vely = groline[52:60]
    velz = groline[60:68]
    return {k: v.strip() for k, v in {
        "resid":resid, 
        "resname":resname, 
        "atomname":atomname, 
        "atomid":atomid, 
        "posx":posx, 
        "posy":posy, 
        "posz":posz, 
        "velx":velx, 
        "vely":vely, 
        "velz":velz}.items()}

class Replacement(object):
    def __init__(self, simulation_box, selection_string, composition_name):
        self.simulation_box = simulation_box
        self.config = self.simulation_box.alchex_config
        self.selection_string = selection_string
        self.composition = self.config.compositions[composition_name]
        self.source_residues = self.simulation_box.mda_universe.select_atoms(selection_string).residues
        self.plan_replacements()
    def plan_replacements(self):
        self.replacements = []
        available = list(range(len(self.source_residues)))
        shuffle(available)
        asc_composition = sorted(self.composition.items(), key=lambda x : x[1])
        for resname, fraction in asc_composition:
            ideal_number = int(round(fraction * len(self.source_residues)))
            res_idxs = []
            if ideal_number >= 1:
                res_idxs = available[:ideal_number]
                available = available[ideal_number:]
            for res_idx in res_idxs:
                self.replacements.append((res_idx, resname))
    def perform_replacements(self):
        self.replace_groups = []
        for res_idx, resname in self.replacements:
            from_structure = self.source_residues[res_idx]
            from_resname = from_structure.atoms[0].resname
            from_residue = ResidueStructure(from_structure, self.config.parameters[from_resname])
            to_residue = self.config.get_reference_structure(resname)
            rep_map = self.config.exchange_maps[from_resname][resname]
            exchange = rep_map.run(from_residue, to_residue, new_resid=from_structure.atoms[0].resid)
            self.replace_groups.append(({"resid": [str(from_structure.atoms[0].resid)]}, exchange))
    def save_grofile(self, output_gro):
        out_lines = []
        rep_made = [False for x in self.replace_groups]
        with open(self.simulation_box.input_filename, "r") as input_fh:
            atoms_started=False
            atom_count = 0
            for line in input_fh:
                writeable = [line.split("\n")[0]]
                d = gro_to_dict(line)
                atoms_started = True
                for idx, (condition, residue) in enumerate(self.replace_groups):
                    added = rep_made[idx]
                    c_sat = True
                    for k, v in condition.items():
                        if (k not in d) or (d[k] not in v):
                            c_sat = False
                            break
                    if c_sat:
                        writeable = []
                        if not rep_made[idx]:
                            writeable = residue.as_gro()
                        rep_made[idx] = True
                        break
                out_lines += writeable
        out_lines[1] = str(len(out_lines) - 3)
        for d in range(2, len(out_lines)-1):
            out_lines[d] = out_lines[d][:15] + str(d-1).rjust(5) + out_lines[d][20:]
        with open(output_gro, "w") as output_fh:
            output_fh.write("\n".join(out_lines))


class SimulationBox(object):
    def __init__(self, alchex_config):
        self.alchex_config = alchex_config
        self.mda_universe = None
        self.replacements = {}
    def load_universe(self, universe_filename):
        self.input_filename = universe_filename
        self.mda_universe = mda.Universe(universe_filename)
    def add_replacement(self, selection_string, composition_name, replacement_name="default"):
        self.replacements[replacement_name] = Replacement(self, selection_string, composition_name)
    def perform_replacement(self, output_file, replacement_name="default"):
        rep = self.replacements[replacement_name]
        rep.perform_replacements()
        rep.apply_modifications_gro(output_file)

'''
popc_residues   = popc_pdb.select_atoms("resname POPC")

lines = []
for popc_residue in popc_residues.residues:
    popc_structure = ResidueStructure(popc_residue, popc)
    swapped = g.run(popc_structure, dlpg_structure)
    lines += swapped.as_pdb()
'''



'''


import MDAnalysis as mda
import sys
sys.path.append("../")
from cgswap.residue import ResidueStructure
from cgswap.gromacs import GromacsITPFile
from cgswap.exchange_map import ExchangeMap

dlpg_itp = GromacsITPFile("../data/DLPG.itp")
dlpg = dlpg_itp.read_residue("DLPG")
dvpe_itp = GromacsITPFile("../data/DVPE.itp")
dvpe = dvpe_itp.read_residue("DVPE")
cdl_itp = GromacsITPFile("../data/CDL.itp")
cdl = cdl_itp.read_residue("CDL0")
chol_itp = GromacsITPFile("../data/chol.itp")
chol = chol_itp.read_residue("CHOL")
dppc_itp = GromacsITPFile("../data/dppc.itp")
dppc = dppc_itp.read_residue("DPPC")
popc_itp = GromacsITPFile("../data/popc_patch/popc_ac.itp")
popc = popc_itp.read_residue("POPC")

dppc_pdb       = mda.Universe("../data/dppc.pdb")
dppc_residue   = dppc_pdb.select_atoms("resname DPP").residues[0]
dppc_structure = ResidueStructure(dppc_residue, dppc)

chol_pdb       = mda.Universe("../data/chol.pdb")
chol_residue   = chol_pdb.select_atoms("resname CHO").residues[0]
chol_structure = ResidueStructure(chol_residue, chol)

dlpg_pdb       = mda.Universe("../data/DLPG-em.gro")
dlpg_residue   = dlpg_pdb.select_atoms("resname DLPG").residues[0]
dlpg_structure = ResidueStructure(dlpg_residue, dlpg)

cdl_pdb       = mda.Universe("../data/CDL0.gro")
cdl_residue   = cdl_pdb.select_atoms("resname CDL").residues[0]
cdl_structure = ResidueStructure(cdl_residue, cdl)

popc_pdb      = mda.Universe("../data/popc_patch/sample.gro")

z = ExchangeMap()
z.new(from_itp=dppc, to_itp=chol, method="martini.static_planar_alignment", draw=False, clusters=[
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
'''
'''
z.run(dppc_structure, chol_structure)

y = ExchangeMap()
y.new(from_itp=dppc, to_itp=dlpg, method="martini.lipid", draw=False)

y.run(dppc_structure, dlpg_structure)

x = ExchangeMap()
x.new(from_itp=popc, to_itp=cdl, method="martini.lipid_to_card", draw=False)

hgs = popc_pdb.select_atoms(x.actions[1]["from_selector"])

from cgswap.geometry import PointCloud, plot_3d

a = PointCloud(3)
a.add_points(hgs.coordinates())
b = PointCloud(3)
b.add_points(hgs.coordinates())
#plot_3d(a)

d = a.find_friends(b, 1)

a_resid1, a_resid2, _ = d[0]

resid1 = hgs.atoms[a_resid1].resid
resid2 = hgs.atoms[a_resid2].resid

popc_residue1   = popc_pdb.select_atoms("resid " + str(resid1))
popc_residue2   = popc_pdb.select_atoms("resid " + str(resid2))
popc_structure1 = ResidueStructure(popc_residue1, popc)
popc_structure2 = ResidueStructure(popc_residue2, popc)

from cgswap.residue import MultiResidueStructure

a = MultiResidueStructure([popc_structure1, popc_structure2])

x.run(a, cdl_structure)
'''
'''
g = ExchangeMap()
g.new(from_itp=popc, to_itp=dlpg, method="martini.lipid", draw=False)

popc_residues   = popc_pdb.select_atoms("resname POPC")

lines = []
for popc_residue in popc_residues.residues:
    popc_structure = ResidueStructure(popc_residue, popc)
    swapped = g.run(popc_structure, dlpg_structure)
    lines += swapped.as_pdb()
with open("../output.pdb", "w") as ofile:
    ofile.write("\n".join(lines))
    '''