from os import path, makedirs
import sys
sys.path.append("/sansom/n15/shil3498/dphil/prj/2016-01-06_Phospholipids/alchex")

from alchex.gromacs_interface import SimulationContainer, GromacsWrapper

gwrap = GromacsWrapper("/sbcb/packages/opt/Linux_x86_64/gromacs/5.1/bin/gmx_avx")

def setup(folder="3ZE3-phosynth"):
    gro_path = path.abspath(path.join(folder, "alchex", "all-em.gro"))
    topol_path = path.abspath(path.join(folder, "alchex", "replaced.top"))
    mdp_path = "/sansom/n15/shil3498/dphil/prj/2016-01-06_Phospholipids/alchex/phospho_synth_sims/replaced/1000ns_run.mdp"
    
    sims = SimulationContainer(path.join(folder,"mdrun"), gwrap)
    protein_atoms = []
    protein_residues = set()
    lipid_atoms = []
    water_sol_atoms = []
    with open(gro_path, "r") as file_handle:
        for line in file_handle:
            if len(line) < 40:
                continue
            resid = line[0:5].strip()
            resname = line[5:10].strip()
            atomid = line[15:20].strip()
            if resname in ["W", "ION"]:
                water_sol_atoms.append(atomid)
            elif resname in ["PI3P","CDDG", "POPE", "PODG", "CDL0", "POPG", "POPA"]:
                lipid_atoms.append(atomid)
            else:
                protein_atoms.append(atomid)
                protein_residues.add(resname)
    with open(sims.resolve_path("index.ndx"), "w") as file_handle:
        file_handle.write("[ Protein ]\n")
        file_handle.write("\n".join(protein_atoms))
        file_handle.write("\n[ Lipid ]\n")
        file_handle.write("\n".join(lipid_atoms))
        file_handle.write("\n[ SOL_ION ]\n")
        file_handle.write("\n".join(water_sol_atoms))
    print folder
    print sims.gromacs.grompp(kwargs={
        "-f" : mdp_path,
        "-c" : gro_path,
        "-p" : topol_path,
        "-n" : "index.ndx",
        "-o" : "mdrun-"+folder[:4]+".tpr",
        "-maxwarn" : 3
        })[0]
    print sims.gromacs.mdrun(kwargs={
        "-deffnm" : "mdrun-"+folder[:4],
        "-nsteps" : 100000000
        })[0]
setup("3ZE3-phosynth")
setup("4WIS-phosynth")
setup("5D91-phosynth")
setup("4Q2E-phosynth")