def is_number(x):
    try:
        float(x)
        return True
    except:
        return False

def martini_atom_similarity(atom1, atom2):
    similarity = 0
    if atom1.attrs["type"] == atom2.attrs["type"]:
        similarity += 1
    else:
        key = ".".join(sorted([atom1.attrs["type"] , atom2.attrs["type"]]))
        differences = {
            "C1.C3" : 0.6,
            "P4.P5" : 0.6
        }
        similarity += differences.get(key, 0)
    if atom1.attrs["atom"] == atom2.attrs["atom"]:
        similarity += 1
    return similarity

ITP_FIELDS = {"atoms":["id",
                      "type",
                      "resnr",
                      "residue",
                      "atom",
                      "cgnr",
                      "charge",
                      "mass"],
            "bonds":["i",
                     "j",
                     "funct",
                     "length",
                     "force"],
            "angles":["i",
                      "j",
                      "k",
                      "funct",
                      "angle",
                      "force"],
            "moleculetype": ["molname",
                             "nrexcl"]
                     }


