/* Header file for C-compatible variables/functions in box.f90 */

/* Set/get the number of boxes (replicas) to simulate */
void box_set_num_boxes(int num_replicas);
void box_get_num_boxes(int num_replicas[1]);

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
/* void box_minimum_image(double vector1[3], double vector[3], double vec_out[3]); */

/* Build link cell data structure from scratch */
