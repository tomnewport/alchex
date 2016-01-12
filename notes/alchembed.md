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

