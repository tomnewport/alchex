from cgswap.residue_parameters import ResidueParameters

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
                file_handle.write(k.ljust(25) + "= " + v








