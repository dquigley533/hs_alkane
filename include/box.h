/* Header file for C-compatible variables/functions in box.f90 */

/* Set/get the number of boxes (replicas) to simulate */
void box_set_num_boxes(int num_replicas);
void box_get_num_boxes(int *num_replicas_out);

/* Initialise/destroy the simulation box(es) */
void box_initialise();
void box_destroy();

/* Set/get the matrix of cell vectors */
void box_set_cell(int ibox, double cell_matrix[3][3]);
void box_get_cell(int ibox, double outmat[3][3]);

/* Convert between absolute and fractional coordinates */
void box_cart_to_frac(int ibox, double vector[3], double vec_out[3]);
void box_frac_to_cart(int ibox, double vector[3], double vec_out[3]);

/* Find the shortest image vector connecting two spheres */
void box_minimum_image(int ibox, double vector1[3], double vector2[3], double vec_out[3]);

/* Build & destroy link cell data structure */
void box_set_link_cell_length(double length);
void box_construct_link_cells(int ibox);
void box_destroy_link_cells();

/* Set flag indicating that use of link cells should be bypassed */
void box_set_bypass_link_cells(int flag);

/* Recalculate the matrix of reciprocal lattice vectors */
void box_update_recipmatrix(int ibox);

/* Calculate volume of simulation cell */
double box_compute_volume(int box);

/* Enforce that box can only change volume isotropically */
void box_set_isotropic(int flag);


/* Set flag indicating that Verlet lists should if link cells bypassed */
void box_set_use_verlet_list(int flag);

/* Set flag indicating if periodic rather open boundaries should be used */
void box_set_pbc(int flag);
