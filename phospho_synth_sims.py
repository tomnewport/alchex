'''
This script adds a centred .gro file to each folder and renames DPP to DPPC
'''

import sys

sys.path.append("../")

from glob import glob
from os import path
from subprocess import check_output
from alchex.replacement import ReplacementSystem

GMX = "/sbcb/packages/opt/Linux_x86_64/gromacs/5.1/bin/gmx_avx"

for pdb_filename in glob("/sansom/n15/shil3498/dphil/prj/2016-01-06_Phospholipids/alchex/phospho_synth_sims/input_data/*/*/coarsegrained-system.pdb"):
    folder, filename = path.split(pdb_filename)
    writefolder = path.join("/sansom/n15/shil3498/dphil/prj/2016-01-06_Phospholipids/alchex/phospho_synth_sims/replaced/", path.split(folder)[1]+"-phosynth")
    check_output("cd "+folder+";"+GMX +" editconf -f "+filename+" -o input.gro -c", shell=True)
    with open(path.join(folder, "input.gro"), "r") as input_fh:
        newlines = input_fh.read().replace("DPP ", "DPPC")
    with open(path.join(folder, "input.gro"), "w") as output_fh:
        output_fh.write(newlines)
    d = ReplacementSystem(
        input_structure_filename=path.join(folder,"input.gro"),
        input_topology_filename=path.join(folder,"topol2.top"),
        root_folder=path.split(folder)[1]+"-phosynth"
    )
    
    d.auto_replace({
    "selection": "resname DPPC and leaflet upper",
    "composition" :{
        "POPA" : 1,
        "PI3P" : 1,
        "CDDG" : 1,
        "PODG" : 1,
        "CDL0" : 1,
        "POPG" : 4,
        "POPE" : 12
        }},
    {
    "selection": "resname DPPC and leaflet lower",
    "composition" :{
        "POPA" : 1,
        "PI3P" : 1,
        "CDDG" : 1,
        "PODG" : 1,
        "CDL0" : 1,
        "POPG" : 4,
        "POPE" : 12
        }},
    )
    