PACKMOL_INDENT = "  "

def indent(string, indent=PACKMOL_INDENT):
    string = indent + string
    return string.replace("\n", "\n"+indent)


'''
From the Packmol manual:

resnumbers 0 

In this case the residue numbers of all residues will 
correspond to the molecule of each type, independently
of the residue numbering of the original pdb file. 
This means that if you pack 10 water molecules and 10
ethanol molecules, the water molecules will be 
numbered 1 to 10, and the ethanol molecules will be 
numbered 1 to 10. 

resnumbers 1 

In this case, the residue numbers of the original pdb 
files are kept unmodified. This means that if you pack
10 proteins of 5 residues, the residue numbers will 
be preserved and, therefore, they will be repeated 
for equivalent residues in each molecule of the same 
protein. 

resnumbers 2 

In this case, the residue numbers of all residues for 
this structure will be numbered sequentially according
to the number of residues that are printed previously
in the same file. This means that if you pack 10 
proteins of 5 residues, there will be residue numbers
ranging from 1 to 50. 

resnumbers 3 

In this case, the numbering of the residues will 
correspond to the sequential numbering of all residues 
in the file. That is, if you pack a protein with 150 
residues followed by 10 water molecules, the water 
molecules will be numbered from 151 to 161. 
'''

RESNUMBER_BYTYPE = 0

RESNUMBER_UNMODIFIED = 1

RESNUMBER_SEQUENTIAL = 2

RESNUMBER_ALLRES_SEQUENTIAL = 3

class WrappedPackmolConstraint(object):
    def __init__(self, constraint_type="", *constraint_args):
        self.constraint_type = constraint_type
        self.constraint_args = constraint_args
    def __str__(self):
        return "{constraint_type} {constraint_args}".format(
                constraint_type = self.constraint_type,
                constraint_args = " ".join(map(str,self.constraint_args))
            )

class WrappedPackmolAtomSelection(object):
    def __init__(self, atoms=[], radius=None, centerofmass=False):
        self.atoms = atoms
        self.constraints = []
        self.radius = radius
    def __str__(self):
        string = "atoms {atom_ids}\n".format(atom_ids = " ".join(map(str,self.atoms)))
        if self.radius:
            string += indent("radius {radius}".format(radius=self.radius)) + "\n"
        string += "\n".join(map(lambda x : indent(str(x)), self.constraints))
        string += "\nend atoms"
        return string
        

class WrappedPackmolStructure(object):
    def __init__(self, filename=None, number=None, resnumbers=None):
        self.filename = filename
        self.number = number
        self.resnumbers = resnumbers
        self.constraints = []
        self.atom_selections = []
    def __str__(self):
        string = "structure {filename}".format(filename=self.filename)+"\n"
        if self.number is not None:
            string += indent("number {number}".format(number=self.number))+"\n"
        if self.resnumbers is not None:
            string += indent("resnumbers {resnumbers}".format(resnumbers=self.resnumbers))+"\n"
        string += "\n".join(map(lambda x : indent(str(x)), self.constraints+self.atom_selections))
        string += "\nend structure"
        return string

class WrappedPackmolInput(object):
    def __init__(self, tolerance=2.0, output="packmol.pdb", filetype="pdb"):
        self.tolerance = tolerance
        self.output = output
        self.filetype = filetype
        self.structures = []
    def __str__(self):
        string = "tolerance {tolerance}".format(tolerance=self.tolerance) + "\n"
        string += "output {output}".format(output=self.output) + "\n"
        string += "filetype {filetype}".format(filetype=self.filetype) + "\n\n"
        string += "\n\n".join(map(str, self.structures))
        return string

from subprocess import Popen, PIPE, STDOUT
from os import path
import MDAnalysis as mda
import numpy
from alchex.geometry import PointCloud, plot_3d, TransformationMatrix, cross_sectional_area_3d
import matplotlib.pyplot as plt

PACKMOL_DIR = "/sansom/n15/shil3498/apps/packmol/packmol"

class WrappedPackmol(object):
    def __init__(self, packmol_executable="packmol"):
        self.packmol_executable = packmol_executable
    def run(self, input, directory="."):
        #packmol_process = Popen(PACKMOL_DIR, stdout=PIPE, stdin=PIPE, stderr=STDOUT, cwd=directory)
        #response = packmol_process.communicate(str(input))
        packmol_process = Popen(PACKMOL_DIR, stdin=PIPE, cwd=directory)
        packmol_process.communicate(str(input))

class VesicleBuilder(object):
    def __init__(self, output="vesicle.pdb", diameter=100, cwd="."):
        self.output = output
        self.diameter = diameter
        self.cwd = cwd
        self.molecules = []
    def add_molecule(self, input_pdb, outer=None, inner=None, center=None, orient_tolerance=0.1):
        universe = mda.Universe(path.join(self.cwd, input_pdb))
        spheres = [(k,universe.select_atoms(v)) for k,v in {-1:inner, 0:center, 1:outer}.items() if v is not None][:2]
        displacement = spheres[0][1].centroid() - spheres[1][1].centroid()
        distance = numpy.linalg.norm(displacement)
        mol_axis = displacement / distance
        points = PointCloud(3)
        points.add_points(universe.select_atoms("resname DPP").coordinates())
        cross_sectional_area_3d(points, mol_axis)
        


vesicle = VesicleBuilder(output="smallvesicle.pdb", diameter=30, cwd="packmol_test")
a = vesicle.add_molecule("molecules/2rlf.pdb", outer="name PO4 and prop z > 10", center="name C4A or name C4B")

'''
thickness = 17
mid = 90
margin = 0.2
top = mid + thickness
bottom = mid - thickness
tmid = mid + margin
bmid = mid - margin

prot_head = WrappedPackmolAtomSelection([374,386,410,422,446,470,482,506])
prot_tail = WrappedPackmolAtomSelection([398,434,458,494,518,530,542,554])

prot_head_const = WrappedPackmolConstraint("outside sphere", 0,0,0,top)
prot_tail_const = WrappedPackmolConstraint("inside sphere", 0,0,0,bottom)

prot_head.constraints.append(prot_head_const)
prot_tail.constraints.append(prot_tail_const)

dppc_head_const = WrappedPackmolConstraint("outside sphere", 0,0,0,top)
dppc_tail_const = WrappedPackmolConstraint("inside sphere", 0,0,0,tmid)

dppc_head = WrappedPackmolAtomSelection([2])
dppc_tail = WrappedPackmolAtomSelection([8,12])

dppc_head.constraints.append(dppc_head_const)
dppc_tail.constraints.append(dppc_tail_const)

dppc_outer = WrappedPackmolStructure("molecules/dppc.pdb", 1500, resnumbers=RESNUMBER_SEQUENTIAL)
dppc_outer.atom_selections.append(dppc_head)
dppc_outer.atom_selections.append(dppc_tail)

dppc_inner_head_const = WrappedPackmolConstraint("inside sphere", 0,0,0, bottom)
dppc_inner_tail_const = WrappedPackmolConstraint("outside sphere", 0,0,0, bmid)

dppc_inner_head = WrappedPackmolAtomSelection([2])
dppc_inner_tail = WrappedPackmolAtomSelection([8,12])

dppc_inner_head.constraints.append(dppc_inner_head_const)
dppc_inner_tail.constraints.append(dppc_inner_tail_const)

dppc_inner = WrappedPackmolStructure("molecules/dppc.pdb", 1000, resnumbers=RESNUMBER_SEQUENTIAL)
dppc_inner.atom_selections.append(dppc_inner_head)
dppc_inner.atom_selections.append(dppc_inner_tail)

e = WrappedPackmolStructure("molecules/2rlf.pdb", 30, resnumbers=RESNUMBER_SEQUENTIAL)
e.atom_selections.append(prot_head)
e.atom_selections.append(prot_tail)

d = WrappedPackmolInput()

d.structures.append(dppc_outer)
d.structures.append(dppc_inner)
d.structures.append(e)


x = WrappedPackmol(PACKMOL_DIR)
x.run(d, "/sansom/n15/shil3498/dphil/prj/2016-01-06_Phospholipids/alchex/packmol_test")
'''