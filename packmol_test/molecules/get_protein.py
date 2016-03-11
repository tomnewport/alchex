import MDAnalysis as mda
import sys


u = mda.Universe(sys.argv[1])

protein = u.select_atoms("not resname W and not resname ION and not resname DPP")

annular = u.select_atoms("resname DPP and same resid as around 10 protein", protein=protein)

W = mda.Writer(sys.argv[2], multiframe=True)
W.write(protein + annular)
W.close()
