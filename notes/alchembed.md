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