# hs_alkane
Fortran2003 code (with C bindings) implementing hard-sphere alkane models as described in the following papers from
the Monson group at the University of Massachusetts.

* [A. P. Malanoski and P. A. Monson, *J. Chem. Phys.* 110, 664 (1999).](https://aip.scitation.org/doi/10.1063/1.478123)
* [M. Cao, A. P. Malanoski, J.W. Schroer, P. A. Monson, *Fluid Ph. Equilibria* 228-229, 75 (2005).](https://pubs.acs.org/doi/abs/10.1021/jp902887w)
* [M. Cao and P. A. Monson, *J. Phys. Chem. B 113*, 13866 (2009).](https://www.sciencedirect.com/science/article/pii/S0378381204003796)

This code was created for the work reported in [Bridgwater & Quigley, "Lattice-switching Monte Carlo method for crystals of flexible 
molecules", *Phys. Rev. E* 90, 063313 (2014)](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.90.063313). These are fun models to
play with so others users may find the code useful/interesting. No guarantee of fitness for any particular purpose is implied by 
that statement.

### Monte Carlo moves

Routines which implement the following Monte Carlo trial moves are included.

* Translation of a randomly selected entire alkane chain - `alkane_translate_chain`
* Rigid rotation of a randomly selected chain about the position of the first bead - `alkane_rotate_chain`
* Rotation around a randomly selected torsion angle along the chain - `alkane_bond_rotate`
* Configurational bias Monte Carlo moves (regrowth of the chain from a randomly selected bead) - `alkane_grow_chain`
* Parrinello-Rahman type changes to the lattice vectors defining the simulation cell - `alkane_box_resize`

### Visualisation

The chain topology and snapshots of bead positions are written to `chain.psf` and `chain.dcd` respectively, suitable for visualisation in [VMD](http://www.ks.uiuc.edu/Research/vmd/). The routines for this in `vis_module.f90` might be useful for others who are stuggling to read/write to these formats in their own code as much as I did!
