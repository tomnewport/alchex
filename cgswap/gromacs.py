from cgswap.residue_parameters import ResidueParameters
import re

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
    def get_includes(self):
        r = []
        for line_idx, line in enumerate(self.lines):
            if "#include" in line:
                r.append((line_idx, line))
        return r
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
       
class GromacsGROFile()




a = GromacsTOPFile()
a.from_file("martini_v2.1.itp")


'''Testing code
a = GromacsMDPFile()

a.from_file("/Users/tom/github_repos/alchex/gromacs_scratch/em.mdp")
a.attrs["test"] = "success"
a.to_file("test.mdp")







'''
