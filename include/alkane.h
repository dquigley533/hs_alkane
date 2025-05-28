/* Header file for C-compatible variables/functions in alkane.f90 */


/* Set, get the number of chains in each simulation box/replica */
void alkane_get_nchains(int *nchains);
void alkane_set_nchains(int nchains);

/* Set, get the number of beads on each chain */
void alkane_get_nbeads(int *nbeads);
void alkane_set_nbeads(int nbeads);

/* Initialise/destroy the alkane module data structures */
void alkane_initialise(void);
void alkane_construct_linked_lists(int ibox);
void alkane_destroy(void);

/* Grow a chain using CBMC */
void alkane_grow_chain(int ichain, int ibox, double *rbfactor, int new_conf, int *ifail);

/* Get a chain */
void alkane_get_chain(int ichain, int ibox, int *nbeads_out,int *d_out, double **rchain_ptr);

/* Translation, rotation of a chain */
void alkane_translate_chain(int ichain, int ibox, double *boltz_out, double vec_out[3]);
void alkane_rotate_chain(int ichain, int ibox, double *boltz_out, double quat[4], int bond);

/* Rotate around a torsion angle */
void alkane_bond_rotate(int ichain, int ibox, double *boltz_out, int *ia, double *angle, int allow_flip);

/* Check if any chain overlaps with some over chain (for consistency checking) */
void alkane_check_chain_overlap(int ibox, int *overlap);

/* Check that a chain has the correct internal geometry (for constency checking) */
void alkane_check_chain_geometry(int ichain, int ibox, int *violated);

/* Update linked lists - for use after accepting any move which changes bead positions */
void alkane_update_linked_lists(int ibead, int ichain, int ibox, double old_pos[3], double new_pos[3]);

/* Construct Verlet lists if we're using those instead of link cells. */
void alkane_construct_neighbour_list(int ibox);

/* Get and set maximum translational displacement */
void alkane_get_dr_max(double *dum_dr);
void alkane_set_dr_max(double dr_max);

/* Get and set maximum rotation angle for rotation moves */
void alkane_get_dt_max(double *dum_dt);
void alkane_set_dt_max(double dt_max);

/* ...as above but for rotations around the first bond on the chain */
void alkane_get_axis_max(double *dum_axis);
void alkane_set_axis_max(double dum_axis);

/* Get and set maximum rotation angle for dihedral moves */
void alkane_get_dh_max(double *dum_dh);
void alkane_set_dh_max(double dh_max);

/* Get and set maximum volume/cell vector change parameter */
void alkane_get_dv_max(double *dum_dv);
void alkane_set_dv_max(double dv_max);

/* Get and set maximum number of trials at each stage of a CBMC regrowth move */
void alkane_get_ktrial(int *dum_ktrial);
void alkane_set_ktrial(int ktrial);

/* Get and set maximum number of beads to regrow in a CBMC move */
void alkane_get_max_regrow(int *dum_max_regrow);
void alkane_set_max_regrow(int dum_max_regrow);

/* Change the simulation box by a matrix delta_H */
void alkane_change_box(int ibox, double delta_H[3][3]);

/* Scale the simulation cell vectors */
void alkane_box_scale(int ibox, double scaleA, double scaleB, double scaleC);

/* Move Carlo move which changes the size of the box */
void alkane_box_resize(double pressure, int ibox, double *boltz_out, int reset);


/* Get and distance between spheres */
void alkane_get_bondlength(double *dumlength);
void alkane_set_bondlength(double L);

void alkane_get_bondangle(double *dumangle);
void alkane_set_bondangle(double bondangle);
