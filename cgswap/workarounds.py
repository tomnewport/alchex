class WAEditableResidue(object):
    def __init__(self, resname, resid):
        self.resname = resname
        self.resid = resid
        self.atoms = []
    def import_mdanalysis_atoms(self, universe, selection_string="all", resids=None):
        print(universe.select_atoms(selection_string))
        print("OK")
