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
        self.movebadrandom = False
        self.discale = 1.1
        self.maxit = 20
    def __str__(self):
        string = "tolerance {tolerance}".format(tolerance=self.tolerance) + "\n"
        if self.movebadrandom:
            string += "movebadrandom\n"
        string += "discale "+str(self.discale)+"\n"
        string += "maxit "+str(self.maxit)+"\n"
        string += "output {output}".format(output=self.output) + "\n"
        string += "filetype {filetype}".format(filetype=self.filetype) + "\n\n"
        string += "\n\n".join(map(str, self.structures))
        return string

from subprocess import Popen, PIPE, STDOUT
from os import path
import MDAnalysis as mda
import numpy
from alchex.geometry import PointCloud, plot_3d, TransformationMatrix, voronoi_shell_3d, cross_sectional_area_3d
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
    def __init__(self, output="vesicle.pdb", radius=100, cwd="."):
        self.output = output
        self.radius = radius
        self.cwd = cwd
        self.molecules = []
    def add_molecule(self, input_pdb, outer=None, inner=None, center=None, orient_tolerance=0.9, ratio=1):
        universe = mda.Universe(path.join(self.cwd, input_pdb))
        spheres = [(k,universe.select_atoms(v)) for k,v in [(-1,inner), (0,center), (1,outer)] if v is not None][:2]
        displacement = spheres[0][1].centroid() - spheres[1][1].centroid()
        distance = numpy.linalg.norm(displacement)
        mol_axis = displacement / distance
        points = PointCloud(3)
        points.add_points(universe.atoms.coordinates())
        area = cross_sectional_area_3d(points, mol_axis)
        self.molecules.append((
            spheres[0],
            spheres[1],
            input_pdb,
            area,
            distance,
            orient_tolerance,
            ratio
        ))
    def generate_input(self):
        input_object = WrappedPackmolInput(output=self.output)  
        outer_unitarea = 0
        inner_unitarea = 0
        for sphere0, sphere1, input_pdb, area, distance, orient_tolerance, ratio in self.molecules:
            if sphere0[0] == -1:
                inner_unitarea += area * ratio
            if sphere1[0] == 1:
                outer_unitarea += area * ratio
        for sphere0, sphere1, input_pdb, area, distance, orient_tolerance, ratio in self.molecules:
            if sphere0[0] == -1 and sphere1[0] == 1:
                sphere_inside_radius = self.radius - ((distance * orient_tolerance)/2)
                sphere_outside_radius = self.radius + ((distance * orient_tolerance)/2)
                unitarea = inner_unitarea
            elif sphere0[0] == -1:
                sphere_inside_radius = self.radius - (distance * orient_tolerance)
                sphere_outside_radius = self.radius
                unitarea = inner_unitarea
            else:
                sphere_inside_radius = self.radius 
                sphere_outside_radius = self.radius + (distance * orient_tolerance)
                unitarea = outer_unitarea
            # unitarea/ is number of units that can be packed in
            # one unit contains <ratio> of this molecule
            total_area = ((4 * 3.141592 * sphere_inside_radius )**2)

            packnumber = (total_area/(unitarea*10)) * ratio
            print packnumber
            #packnumber = 0.008 * total_area
            #packnumber = 6
            structure = WrappedPackmolStructure(filename=input_pdb, number = int(round(packnumber)), resnumbers=3)
            sphere_inside_selection = WrappedPackmolAtomSelection([x.number+1 for x in sphere0[1].atoms])
            sphere_outside_selection = WrappedPackmolAtomSelection([x.number+1 for x in sphere1[1].atoms])
            sphere_inside_constraint = WrappedPackmolConstraint("inside sphere",0,0,0,sphere_inside_radius)
            sphere_outside_constraint = WrappedPackmolConstraint("outside sphere",0,0,0,sphere_outside_radius)
            sphere_inside_selection.constraints.append(sphere_inside_constraint)
            sphere_outside_selection.constraints.append(sphere_outside_constraint)
            structure.atom_selections.append(sphere_inside_selection)
            structure.atom_selections.append(sphere_outside_selection)
            input_object.structures.append(structure)
        input_object.movebadrandom = True
        input_object.discale = 1.5
        input_object.maxit = 40
        print input_object
        return input_object
    def build(self):
        packmol = WrappedPackmol(self.cwd)
        packmol.run(self.generate_input(), self.cwd)









vesicle = VesicleBuilder(output="smallvesicle.pdb", radius=400, cwd="packmol_test")
a = vesicle.add_molecule("molecules/2rlf.pdb", outer="name PO4 and prop z > 10", inner="name PO4 and prop z < 10", orient_tolerance=0.2)
a = vesicle.add_molecule("molecules/dppc.pdb", outer="name PO4", center="name C4A or name C4B",ratio=10)
a = vesicle.add_molecule("molecules/dppc.pdb", inner="name PO4", center="name C4A or name C4B",ratio=10)

for radius in [60]:
    vesicle.radius = radius
    vesicle.output = "vesicle-dev-"+str(radius) + ".pdb"
    vesicle.build()

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
