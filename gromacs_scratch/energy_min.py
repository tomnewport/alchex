grofile = "/sansom/n15/shil3498/dphil/prj/2016-01-06_Phospholipids/phosynth/output.gro"
mdpfile = "/sansom/n15/shil3498/dphil/prj/2016-01-06_Phospholipids/phosynth/gromacs_scratch/em.mdp"
topfile = "/sansom/n15/shil3498/dphil/prj/2016-01-06_Phospholipids/phosynth/gromacs_scratch/topol.top"

grompp -f em.mdp -maxwarn 5 -c output.gro -o temp_em -p sample.top