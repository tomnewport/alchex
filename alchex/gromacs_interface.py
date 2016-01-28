from alchex.residue_parameters import ResidueParameters
import re
import threading
import MDAnalysis as mda
from subprocess import Popen, check_output, CalledProcessError, PIPE, STDOUT
from shutil import copy, copyfile, rmtree, copytree
from os import path, makedirs
from collections import OrderedDict
from copy import deepcopy

'''
class AsyncOSCommand(object):
    def __init__(self, command, cwd="."):
        self.cwd      = cwd
        self.command  = command
        self.thread   = None
        self.waiting  = True
        self.complete = False
    def run(self, callback, callback_args):
        self.waiting = False
        def threaded_run(callback, command, cwd):
            try:
                process = check_output(command, shell=True, cwd=cwd, stderr=STDOUT)
                self.returncode = 0
                self.output = process
            except CalledProcessError as err:
                self.returncode = err.returncode
                self.output = err.output
            callback()
            self.complete = True
            return
        def _callback():
            callback(*callback_args)
        self.thread = threading.Thread(target = threaded_run, args=(callback, self.command, self.cwd))
        self.thread.start()
'''
class GromacsWrapper(object):
    def __init__(self, gromacs_executable):
        self._exec_gmx = gromacs_executable
        self._cwd = "."
    def run(self):
        batch = []
        for idx, command in enumerate(self._todo):
            if command.waiting:
                command.run()
    def cd(self, cwd):
        self._cwd = path.abspath(path.join(self._cwd, cwd))
    def __getattr__(self, name):
        def unknown_call(args=[], kwargs={}):
            return self.generic_gromacs_call(name, args, kwargs)
        return unknown_call
    def arg_encode(self, args):
        argstring = ""
        for argname, argvalue in args.items():
            if argvalue is not None and argvalue is not True:
                argstring += "{name} {value} ".format(name=argname, value=argvalue)
            elif argvalue is True:
                argstring += "{name} ".format(name=argname)
        return argstring[:-1]
    def _stdin(self,stdin):
        if stdin != []:
            return 'echo "' + " ".join([str(x) for x in stdin]) + '" | '
        else:
            return ""
    def shell(self, process, args=[], kwargs={}):
        command = self.build_command(process, args=args, kwargs=kwargs)
        try:
            output = check_output(command, shell=True, cwd=self._cwd, stderr=STDOUT)
            retcode = 0
        except CalledProcessError as cpe:
            output = cpe.output
            retcode = cpe.returncode
        return retcode, output
    def build_command(self, process, args=[], kwargs={}):
        stdin = ""
        if "_stdin" in kwargs:
            stdin = self._stdin(kwargs["_stdin"])
        command = "{stdin}{process} {args} {kwargs}".format(
            stdin = stdin, 
            process = process,
            args=" ".join([str(x) for x in args]), 
            kwargs=self.arg_encode(kwargs)
            )
        return command
    def generic_gromacs_call(self, command, args, kwargs):
        args = [command] + list(args)
        return self.shell(self._exec_gmx, args=args, kwargs=kwargs)

class SimulationContainer(object):
    def __init__(self, root_path, gromacs_wrapper):
        self._root = path.abspath(root_path)
        self.gromacs = gromacs_wrapper
        self.makedirs("/")
        self.cd()
        self.universes = {}
    def universe(self, filename):
        if filename not in self.universes:
            self.universes[filename] = mda.Universe(self.resolve_path(filename))
        return self.universes[filename]
    def delete_folder(self, delpath="/"):
        rmtree(self.resolve_path(delpath))
    def add_file(self, source_file, destination="/", rename=False):
        destination_folder = self.makedirs(destination)
        if rename is False:
            copy(source_file, destination_folder)
        else:
            copyfile(source_file, path.join(destination_folder, rename))
    def copy_folder(self, source_path, destination_path):
        source_path      = self.resolve_path(source_path)
        destination_path = self.resolve_path(destination_path)
        copytree(source_path, destination_path)
    def copy_file(self, source_path, destination_path, rename=False):
        source_path      = self.resolve_path(source_path)
        destination_path = self.resolve_path(destination_path)
        if rename is False:
            copy(source_path, destination_path)
        else:
            copyfile(source_path, path.join(destination_path, rename))
    def makedirs(self, dirpath):
        path_to_make = self.resolve_path(dirpath)
        if not path.exists(path_to_make):
            makedirs(path_to_make)
        return path_to_make
    def cd(self, newpath="/"):
        self.gromacs.cd(self.resolve_path(newpath))
    def resolve_path(self, newpath):
        if newpath == "/":
            return self._root
        elif newpath[0] == "/":
            return path.join(self._root, newpath[1:])
        else:
            return path.join(self.gromacs._cwd, newpath)

'''
c = SimulationContainer("testcontainer", GromacsWrapper(_exec_gmx="/sbcb/packages/opt/Linux_x86_64/gromacs/5.1/bin/gmx_sse"))
c.delete_folder()
c.add_file("/sansom/n15/shil3498/dphil/prj/2016-01-06_Phospholipids/phosynth/gromacs_scratch/single_dlpg/dlpg.gro")
c.add_file("/sansom/n15/shil3498/dphil/prj/2016-01-06_Phospholipids/phosynth/gromacs_scratch/single_dlpg/topol.top")
c.add_file("/sansom/n15/shil3498/dphil/prj/2016-01-06_Phospholipids/phosynth/gromacs_scratch/single_dlpg/martini_v2.1.itp")
c.add_file("/sansom/n15/shil3498/dphil/prj/2016-01-06_Phospholipids/phosynth/gromacs_scratch/single_dlpg/em.mdp")
c.makedirs("em")
c.copy_file("em.mdp","em")
c.copy_file("dlpg.gro","em")
c.copy_file("martini_v2.1.itp", "em")
c.copy_file("topol.top", "em")
c.copy_file("em.mdp", "em")
c.makedirs("em/run")
c.cd("em")
c.gromacs.grompp(kwargs={"-f":"em.mdp", "-c":"dlpg.gro", "-p":"topol.top", "-o":"run/em.tpr","-maxwarn":"1"})
c.cd("run")
c.gromacs.mdrun(kwargs={"-deffnm":"em"})
'''
'''
a = GromacsWrapper(_exec_gmx="/sbcb/packages/opt/Linux_x86_64/gromacs/5.1/bin/gmx_sse")
a.cd("/sansom/n15/shil3498/dphil/prj/2016-01-06_Phospholipids/phosynth/gromacs_scratch/single_dlpg")
a.shell("mkdir test")
a.grompp(kwargs={"-f":"em.mdp", "-c":"dlpg.gro", "-p":"topol.top", "-o":"test/em.tpr","-maxwarn":"1"})[1]
a.cd("test")
a.mdrun(kwargs={"-deffnm":"em"})
'''


class GromacsITPFile(object):
    def __init__(self, filename):
        self.filename = filename
    def read_residue(self, residue_name):
        colnames = {"atoms" : "id type resnr residue atom cgnr charge".split(),
                   "bonds": "i  j   funct   length  force.c.".split(),
                   "angles": "i  j  k   funct   angle   force.c.".split(),
                   "molname": "molname       nrexcl".split()}
        tables = {}
        with open(self.filename, "r") as file_handle:
            sectionline = 0
            in_table = False
            correct_section = False
            for line in file_handle:
                stripline = line.strip()
                if len(stripline) > 0:
                    if stripline[0] == "[" and stripline[-1] == "]":
                        table_name = stripline[1:-1].strip()
                        sectionline = 0
                        in_table = True
                        if table_name in colnames:
                            table_columns = colnames[table_name]
                        else:
                            table_columns = None
                    elif stripline[0] == ";" and sectionline == 1:
                        pass
                        #table_columns = line[1:].split()
                    elif in_table:
                        if table_name == "moleculetype":
                            if stripline.split()[0] == residue_name:
                                correct_section = True
                            else:
                                correct_section = False
                        if table_name in ["atoms", "bonds", "angles"] and correct_section and stripline[0] not in ["#", ";"]:
                            if table_name not in tables:
                                tables[table_name] = []
                            tables[table_name].append({key: value for key, value in zip(table_columns, stripline.split())})
                        
                sectionline += 1
        rp = ResidueParameters(residue_name)
        for atom in tables["atoms"]:
            rp.add_atom(atom["id"], atom)
        for bond in tables["bonds"]:
            rp.add_bond(bond["i"], bond["j"], bond)
        for angle in tables["angles"]:
            rp.add_angle(angle["i"], angle["j"], angle["k"], angle)
        return rp

class GromacsMDPFile(object):
    def __init__(self):
        self.attrs = {}
    def from_file(self, filename):
        with open(filename, "r") as file_handle:
            for line in file_handle:
                if "=" in line and line[0] != ";":
                    k, v = line.split("=")
                    self.attrs[k.strip()] = v.strip()
    def to_file(self, filename):
        with open(filename, "w") as file_handle:
            for k, v in self.attrs.items():
                file_handle.write(k.ljust(25) + "= " + str(v) + "\n")

class GromacsTOPFile(object):
    def __init__(self):
        self.lines = []
    def from_file(self, filename):
        with open(filename) as file_handle:
            file_data = file_handle.read()
        self.lines = file_data.split("\n")
    def to_file(self, filename):
        with open(filename, "w") as file_handle:
            file_handle.write("\n".join(self.lines))
    def get_includes(self):
        r = []
        for line_idx, line in enumerate(self.lines):
            if "#include" in line:
                r.append((line_idx, line, re.findall('"([^"]*)"', line)[0]))
        return r
    def modify_molecules(self, new_molecules):
        table_head = re.compile(r'^\s*\[\s*molecules\s*\]\s*$')
        startline = None
        endline = None
        for line_id, line in enumerate(self.lines):
            if startline is None:
                if table_head.search(line):
                    startline = line_id
            else:
                stripline = line.strip()
                if len(stripline) == 0 or stripline[0] == "[":
                    endline = line_id
        newlines = [resname.ljust(5) + str(rescount) for resname, rescount in new_molecules]
        self.lines = self.lines[:startline+1] + newlines + self.lines[endline:]
    def get_tables(self):
        table_head = re.compile(r'^\s*\[\s*(\w[\w\s]\w*)\s*\]\s*$')
        tables = []
        table_name = None
        for line_id, line in enumerate(self.lines):
            th_search = table_head.search(line)
            if th_search is not None:
                table_name = th_search.group(1)
                if len(tables) > 0:
                    tables[-1][1] = line_id
                tables.append([line_id, len(self.lines)-1, table_name, [], []])
            elif table_name is not None:
                tables[-1][3].append((line_id, line))
                stripline = line.strip()
                if len(stripline) > 0:
                    if stripline[0] not in [";","#"]:
                        tables[-1][4].append(stripline.split())
        return tables
    def get_residue(self, target_resname):
        tables = self.get_tables()
        rtables = []
        for table in tables:
            start_line, end_line, table_name, table_raw, table_parsed = table
            if table_name == "moleculetype":
                resname = table_parsed[0][0]
                if resname == target_resname:
                    rtables.append(table)
                elif len(rtables) != 0:
                    return rtables
            elif len(rtables) > 0:
                rtables.append(table)
    def residue_parameters(self, target_resname):
        colnames = {"atoms" : "id type resnr residue atom cgnr charge".split(),
                   "bonds": "i  j   funct   length  force.c.".split(),
                   "angles": "i  j  k   funct   angle   force.c.".split(),
                   "molname": "molname       nrexcl".split()}

        rp = ResidueParameters(target_resname)
        restable = self.get_residue(target_resname)


        for start_line, end_line, table_name, table_raw, table_parsed in restable:
            if table_name == "atoms":
                rows = [
                        {k:v for k, v in zip(colnames["atoms"], row)} 
                    for row in table_parsed if len(row) == len(colnames["atoms"])
                    ]
                for atom in rows:
                    rp.add_atom(atom["id"], atom)
            elif table_name == "bonds":
                rows = [
                        {k:v for k, v in zip(colnames["bonds"], row)} 
                    for row in table_parsed if len(row) == len(colnames["bonds"])
                    ]
                for bond in rows:
                    rp.add_bond(bond["i"], bond["j"], bond)
        return rp


TOPFILE_FIELDS = {
    "moleculetype" : ["name", "nrexcl"],
    "atoms"        : ["id", "type", "resnr", "residue", "atom", "cgnr", "charge", "mass"],
    "bonds"        : ["ai", "aj", "funct", "c0", "c1"],
    "angles"       : ["ai", "aj", "ak", "funct", "c0", "c1"],
    "dihedrals"    : ["ai", "aj", "ak", "al", "funct", "c0", "c1", "c2"],
    "pairs"        : ["ai", "aj", "funct", "c1", "c1"],
    "system"       : ["system_name"],
    "molecules"    : ["moltype", "count"]
}

MOLTYPE_TABLES = {
    "atoms", "bonds", "angles", "dihedrals", "pairs"
}

def gromacs_top_table_string(name, rows):
    if len(rows) == 0:
        return ""
    else:
        string = "[ {name} ]\n".format(name=name)
        for row in rows:
            string += "\t".join(row.values()) + "\n"
        return string + "\n"
        

class GromacsTOPFileMoltype(object):
    def __init__(self):
        self.atoms = []
        self.name = ""
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.pairs = []
        self.moleculetype = []
    def renumber_atoms(self, renumber_map):
        for row in self.atoms:
            row["id"] = renumber_map[row["id"]]
        for row in self.bonds:
            row["ai"] = renumber_map[row["ai"]]
            row["aj"] = renumber_map[row["aj"]]
        for row in self.angles:
            row["ai"] = renumber_map[row["ai"]]
            row["aj"] = renumber_map[row["aj"]]
            row["ak"] = renumber_map[row["ak"]]
        for row in self.dihedrals:
            row["ai"] = renumber_map[row["ai"]]
            row["aj"] = renumber_map[row["aj"]]
            row["ak"] = renumber_map[row["ak"]]
            row["al"] = renumber_map[row["al"]]
        for row in self.pairs:
            row["ai"] = renumber_map[row["ai"]]
            row["aj"] = renumber_map[row["aj"]]
    def renumber_charge_groups(self, renumber_map):
        for row in self.atoms:
            row["cgnr"] = renumber_map[row["cgnr"]]
    def renumber_residues(self, renumber_map):
        for row in self.atoms:
            row["resnr"] = renumber_map[row["resnr"]]
    def topcat(self, cattop):
        cattop = deepcopy(cattop)
        if len(self.atoms) > 0:
            last_atom_id = self.atoms[-1]["id"]
            last_residue_id = self.atoms[-1]["resnr"]
            last_charge_group = self.atoms[-1]["cgnr"]
        else:
            last_atom_id = 0
            last_residue_id = 0
            last_charge_group = 0
        cattop_atom_ids = [x["id"] for x in cattop.atoms]
        cattop_res_ids = [x["resnr"] for x in cattop.atoms]
        cattop_charge_groups = [x["cgnr"] for x in cattop.atoms]
        cattop.renumber_atoms({x: str(int(x) + int(last_atom_id)) for x in cattop_atom_ids})
        cattop.renumber_residues({x: str(int(x) + int(last_residue_id)) for x in cattop_res_ids})
        cattop.renumber_charge_groups({x: str(int(x) + int(last_charge_group)) for x in cattop_charge_groups})
        self.atoms += cattop.atoms
        self.bonds += cattop.bonds
        self.angles += cattop.angles
        self.dihedrals += cattop.dihedrals
        self.pairs += cattop.pairs
    def as_file(self):
        return (
              gromacs_top_table_string("moltype", [self.moleculetype])
            + gromacs_top_table_string("atoms", self.atoms)
            + gromacs_top_table_string("bonds", self.bonds)
            + gromacs_top_table_string("angles", self.angles)
            + gromacs_top_table_string("dihedrals", self.dihedrals)
            + gromacs_top_table_string("pairs", self.pairs)
            )
            

class GromacsEditableTOPFile(object):
    '''
    This TOPFile class can read and modify
    gromacs topologies. It should eventually
    be merged with the main TOPFile class but
    needs:
    get_includes       - list included scripts
    residue_parameters - return a residue_parameters
                         object
    modify_molecules   - modify molecules table
    '''
    def __init__(self):
        self.tables = OrderedDict()
        self.moltypes = OrderedDict()
    def preprocess_line(self, line):
        if ";" in line:
            line = line.split(";")[0]
        line = line.strip()
        return line
    def from_file(self, filename):
        table_header = re.compile('\[\s*(\w*)\s*\]')
        with open(filename, "r") as file_handle:
            for line in file_handle:
                line = self.preprocess_line(line)
                if line != "":
                    th = table_header.search(line)
                    if th is not None:
                        table_name = th.group(1)
                        if table_name in TOPFILE_FIELDS:
                            table_columns = TOPFILE_FIELDS[table_name]
                        else:
                            table_columns = None
                        if table_name == "moleculetype":
                            moltype = GromacsTOPFileMoltype()
                        elif table_columns is not None and not table_name in MOLTYPE_TABLES:
                            self.tables[table_name] = []                            
                    elif table_columns is not None:
                        table_cells = line.split()
                        table_cells += ["","","","",""]
                        table_cells = table_cells[:len(table_columns)]
                        table_row = OrderedDict(zip(table_columns, table_cells))
                        if table_name == "moleculetype":
                            moltype.moleculetype = table_row
                            moltype.name = table_cells[0]
                            self.moltypes[moltype.name] = moltype
                        elif table_name == "atoms":
                            moltype.atoms.append(table_row)
                        elif table_name == "bonds":
                            moltype.bonds.append(table_row)
                        elif table_name == "angles":
                            moltype.angles.append(table_row)
                        elif table_name == "dihedrals":
                            moltype.dihedrals.append(table_row)
                        elif table_name == "pairs":
                            moltype.pairs.append(table_row)
                        else:
                            self.tables[table_name].append(table_row)
    def to_file(self, filename):
        with open(filename, "w") as file_handle:
            for moltype in self.moltypes.values():
                file_handle.write(moltype.as_file())
            for k, v in self.tables.items():
                file_handle.write(gromacs_top_table_string(k, v))
    def unify_moltype(self, name="unified", nrexcl="1"):
        unified = GromacsTOPFileMoltype()
        unified.moleculetype = OrderedDict([("name", name), ("nrexcl", nrexcl)])
        unified.name = name
        for line in self.tables["molecules"]:
            moltype_name = line["moltype"]
            moltype_count = int(line["count"])
            for _ in range(moltype_count):
                unified.topcat(self.moltypes[moltype_name])
        return unified
      