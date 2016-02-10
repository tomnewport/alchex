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
            string += indent("resnumbers {resnumbers}".format(resnumbers=self.resnumber))+"\n"
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

PACKMOL_DIR = "/sansom/n15/shil3498/apps/packmol/packmol"

class WrappedPackmol(object):
    def __init__(self, packmol_executable="packmol"):
        self.packmol_executable = packmol_executable
    def run(self, input):
        packmol_process = Popen(PACKMOL_DIR, stdout=PIPE, stdin=PIPE, stderr=STDOUT)
        response = packmol_process.communicate(str(input))
        print response[0]

a = WrappedPackmolConstraint("fixed", 0,0,0,21,39,39)
b = WrappedPackmolAtomSelection([1,2,3],1)
c = WrappedPackmolStructure("test.pdb", 10)
d = WrappedPackmolInput()
b.constraints.append(a)
c.atom_selections.append(b)
c.constraints.append(a)
d.structures.append(c)
c.atom_selections.append(b)
c.constraints.append(a)
d.structures.append(c)


x = WrappedPackmol(PACKMOL_DIR)
x.run(d)