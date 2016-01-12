Proposal: Use Jeffreys et Al.'s method (Alchembed) to change lipids.

### Chat with PF

Delete old lipids first, then replace with new lipids. Use Alchembed to resolve clashes.

So... how does alchembed work?
------------------------------

Let's look at `try-alchembed.sh`. 

```bash
# First, prepare a TPR file for energy minimisation
grompp   -f common-files/em-$ff.mdp\
         -c common-files/$protein-$ff.pdb\
         -p common-files/$protein-$ff.top\
         -n common-files/$protein-$ff.ndx\
         -po $protein/$ff/$protein-$ff-em.mdp\
         -o $protein/$ff/$protein-$ff-em\
         -maxwarn 1
```

So we need:

- **-f** Parameters (as mdp)
- **-c** A structure (as pdb)
- **-p** A topology (as top)
- **-n** An index (as ndx)

And we get:

- **-po** An output mdp
- **-o** An output tpr

The energy minimisation then gets run:

```bash
mdrun_d  -deffnm $protein/$ff/$protein-$ff-em\
         -ntmpi 1
```

It then gets re-grompped:

```bash
# Now, prepare the ALCHEMBED TPR file
grompp   -f common-files/alchembed-$ff.mdp\
         -c $protein/$ff/$protein-$ff-em.gro\
         -p common-files/$protein-$ff.top\
         -n common-files/$protein-$ff.ndx\
         -po $protein/$ff/$protein-$ff-alchembed.mdp\
         -o $protein/$ff/$protein-$ff-alchembed\
         -maxwarn 2
```

And run:

```bash
# ..and run on a single core. 
mdrun    -deffnm $protein/$ff/$protein-$ff-alchembed\
         -ntmpi 1
```

The important parts are all within MDP files.

Let's look at em-cg.mdp:

```python
integrator               = steep
emstep                   = 0.001
emtol                    = 1000
nsteps                   = 100
```

Now let's look at alchembed-cg.mdp:

```python
# General simulation stuff
integrator           = md           
define               = -DPOSRES
tinit                = 0.0          
dt                   = 0.01
nsteps               = 1000
nstxout              = 100000
nstvout              = 100000
nstfout              = 10
nstlog               = 10
nstenergy            = 10000
nstxtcout            = 10
xtc_precision        = 1000
coulombtype          = Reaction_Field
rlist                = 1.2
rcoulomb_switch      = 0.0  
rcoulomb             = 1.2
epsilon_r            = 15
epsilon_rf           = 0   
vdw_type             = cutoff 
rvdw_switch          = 0.9
rvdw                 = 1.2
cutoff-scheme        = verlet
coulomb-modifier     = potential-shift-verlet
vdw-modifier         = potential-shift-verlet
verlet-buffer-drift  = 0.005

# Temperature coupling
tcoupl               = Berendsen
tc-grps              = SYSTEM
tau_t                = 1.0 
ref_t                = 310 
gen_vel              = yes
gen_temp             = 310
gen_seed             = 111

# Free energy (alchembed part)
free_energy          = yes
init_lambda          = 0.00
delta_lambda         = 1e-3
sc-alpha             = 0.1000
sc-power             = 1
sc-r-power           = 6
couple-moltype       = PROTEIN
couple-lambda0       = none
couple-lambda1       = vdw
```

Important parts seem to be the free energy:

`free_energy          = yes`

`init_lambda          = 0.00`

`delta_lambda         = 1e-3`

`sc-alpha             = 0.1000`

`sc-power             = 1`

`sc-r-power           = 6`

`couple-moltype       = PROTEIN`

`couple-lambda0       = none`

`couple-lambda1       = vdw`