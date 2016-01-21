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

z.run(dppc_structure, chol_structure,  plot=False)

y = ExchangeMap()
y.new(from_itp=dppc, to_itp=dlpg, method="martini.lipid", draw=False)

a = y.run(dppc_structure, dlpg_structure, plot=False)


x = ExchangeMap()
x.new(from_itp=popc, to_itp=cdl, method="martini.lipid_to_card", draw=True)

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
print(popc_structure1, popc_structure2)
a = MultiResidueStructure([popc_structure1, popc_structure2])

cdl = x.run(a, cdl_structure, plot=True)

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

print("\n".join(cdl.as_gro()))
