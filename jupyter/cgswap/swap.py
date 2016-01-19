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

dppc_pdb       = mda.Universe("../data/dppc.pdb")
dppc_residue   = dppc_pdb.select_atoms("resname DPP").residues[0]
dppc_structure = ResidueStructure(dppc_residue, dppc)

chol_pdb       = mda.Universe("../data/chol.pdb")
chol_residue   = chol_pdb.select_atoms("resname CHO").residues[0]
chol_structure = ResidueStructure(chol_residue, chol)

z = ExchangeMap()
z.new(from_itp=dppc, to_itp=chol, method="martini.static_planar_alignment", draw=False, clusters=[
        [
            ["PO4", "NC3", "GL1"],
            ["ROH"],
            1
        ],
        [
            ["C3A"],
            ["R2", "C2"],
            1
        ],
        [
            ["C3B"],
            ["R3","C2"],
            1
        ]
    ])

z.run(dppc_structure, chol_structure)