/* Header file for C-compatible variables/functions in vis_module.f90 */


/* Write a psf file (topology) */
void write_psf(int nbeads, int nchains);

/* Write charmm-style dcd files */
void write_dcd_header(int nbeads, int nchains);
void write_dcd_snapshot();

