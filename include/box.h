/* Header file for C-compatible variables/functions in box.f90 */


/* Set/get the number of boxes (replicas) to simulate */
void box_set_num_boxes(int num_replicas);
void box_get_num_boxes(int num_replicas[1]);

/* Initialise the simulation box(es) */
void box_initialise();


/* Set/get the matrix of cell vectors */
void box_set_cell(int ibox, double cell_matrix[3][3]);
void box_get_cell(int ibox, double outmat[3][3]);
