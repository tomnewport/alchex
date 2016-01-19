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

From the gromacs manual:

`free_energy          = yes`

Interpolate between topology A (lambda=0) to topology B (lambda=1) and write the derivative of the Hamiltonian with respect to lambda (as specified with dhdl-derivatives), or the Hamiltonian differences with respect to other lambda values (as specified with foreign-lambda) to the energy file and/or to dhdl.xvg, where they can be processed by, for example g_bar. The potentials, bond-lengths and angles are interpolated linearly as described in the manual. When sc-alpha is larger than zero, soft-core potentials are used for the LJ and Coulomb interactions.

`init_lambda          = 0.00`

Starting value for lambda (float). Generally, this should only be used with slow growth (i.e. nonzero delta-lambda). In other cases, init-lambda-state should be specified instead. Must be greater than or equal to 0.

`delta_lambda         = 1e-3`

increment per time step for lambda

`sc-alpha             = 0.1000`

The soft-core alpha parameter, a value of 0 results in linear interpolation of the LJ and Coulomb interactions

`sc-power             = 1`

The power for lambda in the soft-core function, only the values 1 and 2 are supported

`sc-r-power           = 6`

The power of the radial term in the soft-core equation. Possible values are 6 and 48. 6 is more standard, and is the default. When 48 is used, then sc-alpha should generally be much lower (between 0.001 and 0.003).

`couple-moltype       = PROTEIN`

Here one can supply a molecule type (as defined in the topology) for calculating solvation or coupling free energies. There is a special option system that couples all molecule types in the system. This can be useful for equilibrating a system starting from (nearly) random coordinates. free-energy has to be turned on. The Van der Waals interactions and/or charges in this molecule type can be turned on or off between lambda=0 and lambda=1, depending on the settings of couple-lambda0 and couple-lambda1. If you want to decouple one of several copies of a molecule, you need to copy and rename the molecule definition in the topology.

`couple-lambda0       = none`

The Van der Waals interactions are turned off and the charges are zero at lambda=0; soft-core interactions will be required to avoid singularities.

`couple-lambda1       = vdw`

The charges are zero (no Coulomb interactions) at lambda=0

**So in conclusion.**

It looks like our major parameter is couple-moltype. This is defined in the topology, which could make things interesting.

The Problem
===========

I have a population of molecules which I want to replace with a second population of molecules.

Limitations
-----------

- Alchembed must fade in an entire molecule type as defined in the topology file, Therefore at least one alchembed operation is required per residue type.

Solution
--------

1. Source molecules are mapped to destination molecules.
1. For each residue type, starting with the largest:
    1. The residues to be faded in will be designated resname ALCH
    1. While there are still residues outside cutoff distance of selected residues:
        1. A source residue group will be selected.
        1. The ALCH residue is aligned to the source residue
        1. The source residue is deleted
    1. If residues are still available a second residue type may be added.
    1. Alchembed is then run and the process is repeated until all residues are converted.

The most basic grammar required is as follows:

```json
{
"replacements":
    {"POPC" : {
        "CARD" : .2,
        "DPPG" : .1
        }
    },
"verify":{
    "CARD" : .2,
    "DPPG" : .1,
    "POPC" : .7
    }
}
```

More things to consider are:

- More advanced mappings, including
    - n:n mapping
    - alignment

So an improved grammar would be:

``json
{
    "replacements" : 
    [
        {
            "source" : 
                {
                    "selector" : "resname DPPC",
                    "ratio"    : 1
                },
            "destination" : 
                {
                    "selector" : "resname CARD",
                    "ratio"    : 2
                },
            "alignment"  : 
                {
                    "type"     : ("volume" or "regression")
                    "point1"   : {
                        "source": "name GLH",
                        "destination" : "name GLH"
                    }
                    (etc.)
                }
        }
    ]
}
```

For lipids this could be made into a shorthand with ease.

As a worst-case fallback, several isolated points may need to be deleted. This is not ideal, but it isn't a disaster either.