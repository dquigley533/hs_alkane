/* -*- mode: C */

%define DOCSTRING
"Fortran2003 code (with C/Python bindings) implementing hard-sphere alkane models as described
 in the following papers from the Monson group at the University of Massachusetts.

 A. P. Malanoski and P. A. Monson, J. Chem. Phys. 110, 664 (1999).
 M. Cao, A. P. Malanoski, J.W. Schroer, P. A. Monson, Fluid Ph. Equilibria 228-229, 75 (2005).
 M. Cao and P. A. Monson, J. Phys. Chem. B 113, 13866 (2009)."
%enddef

/* alkane.i */
/* N.B. Implementing module docstring using the method described 
at http://swig.org/Doc1.3/Python.html#Python_nn66 stops distutils
from recognising the module name.... 
%module(docstrig=DOCSTRING) alkane
*/
%module alkane
%{
#define SWIG_FILE_WITH_INIT

/* These will be included in the wrapper code */
#include "timer.h"
#include "random.h"
#include "quaternion.h"
#include "box.h"
#include "alkane.h"
%}

/* Standard typemaps */
%include "typemaps.i"

/* Numpy array typemap */
%include "numpy.i"


%init %{
  import_array();
%}

/* Map array array onto arguments used in the C interface */

/* Vectors */
%apply(double IN_ARRAY1[ANY]) {(double vector1[3])};
%apply(double IN_ARRAY1[ANY]) {(double vector2[3])};
%apply(double IN_ARRAY1[ANY]) {(double unit_vector1[3])};
%apply(double IN_ARRAY1[ANY]) {(double unit_vector2[3])};
%apply(double IN_ARRAY1[ANY]) {(double rotation_axis[3])};
%apply(double IN_ARRAY1[ANY]) {(double vector[3])};
%apply(double ARGOUT_ARRAY1[ANY]) {(double vec_out[3])};

%apply(double IN_ARRAY1[ANY]) {(double old_pos[3])};
%apply(double IN_ARRAY1[ANY]) {(double new_pos[3])};


/* Quaternions */
%apply(double ARGOUT_ARRAY1[ANY]){(double quat[4])};
%apply(double IN_ARRAY1[ANY]){(double quaternion[4])};
%apply(double IN_ARRAY1[ANY]){(double quaternion1[4])};
%apply(double IN_ARRAY1[ANY]){(double quaternion2[4])};
%apply(double ARGOUT_ARRAY1[ANY]){(double quat_product[4])};
%apply(double ARGOUT_ARRAY1[ANY]){(double quat_inverse[4])};

/* integers */
%apply int *OUTPUT {int *num_replicas_out};
%apply int *OUTPUT {int *nchains};
%apply int *OUTPUT {int *nbeads};
%apply int *OUTPUT {int *ifail};
%apply int *OUTPUT {int *ia};
%apply int *OUTPUT {int *overlap};
%apply int *OUTPUT {int *violated};

/* doubles */
%apply double *OUTPUT {double *rbfactor};
%apply double *OUTPUT {double *boltz_out};
%apply double *OUTPUT {double *angle};

/* Matrices of cell vectors */
%apply(double IN_ARRAY2[ANY][ANY]) {(double cell_matrix[3][3])};
/*%apply(double ARGOUT_ARRAY2[ANY][ANY]) {(double outmat[3][3])};*/
%apply(int* DIM1, int* DIM2, double** ARGOUTVIEW_ARRAY2) {(int *d1, int *d2, double **outmat_ptr)};



/* Chains */
%apply(int* DIM1, int* DIM2, double** ARGOUTVIEW_ARRAY2) {(int *nbeads_out,int *d_out, double **rchain_ptr)};




/*%apply(double ARGOUT_ARRAY2[ANY][ANY]) {double rchain[3][]}; */

//%apply( int DIM1, double* IN_ARRAY1, int DIM1, double* IN_ARRAY1, int DIM1, double* ARGOUT_ARRAY1 ) {(int d1, double *v1, int d2, double *v2, int dq, double *quat)};

/* Docstring information for compute_neighbour_list */
%feature("autodoc", "timer_init()") timer_init;
%define tmr_ini_str
"
    Initialises the internal calculation timer to zero. This should
    be invoked at the start of any calculation using this module.

    The timer is used to report on performance (MC steps / second)
    and to shutdown cleanly within the walltime of a batch processed
    calculation.

    See also timer_elapsed_time() and timer_check_runtime(). 

    Parameters
    ----------

    (none)

"
%enddef
%feature("docstring", tmr_ini_str) timer_init;

%feature("autodoc", "timer_elapsed_time()") timer_elapsed_time;
%define tmr_elp_str
"
    Returns a double precision number indicating the elapsed
    seconds since timer_init() was invoked.

    See also timer_init() and timer_check_runtime(). 

    Parameters
    ----------

    (none)

    Returns
    -------

    Time in seconds.

"
%enddef
%feature("docstring", tmr_elp_str) timer_elapsed_time;

%feature("autodoc", "timer_elapsed_time()") timer_elapsed_time;
%define tmr_elp_str
"
    Returns a double precision number indicating the elapsed
    seconds since timer_init() was invoked.

    See also timer_init() and timer_check_runtime(). 

    Parameters
    ----------

    (none)

    Returns
    -------

    Time in seconds.

"
%enddef
%feature("docstring", tmr_elp_str) timer_elapsed_time;

%feature("autodoc", "timer_check_continuation()") timer_check_continuation;
%define tmr_chk_str
"
    Returns an integer flag (1/0) indicating if calculation
    is safe to continue without running out of walltime.

    See also timer_init() and timer_elapsed_time(). 

    Returns
    ----------

    1 if remaining walltime > timer_closetime
    0 otherwise (calculation should halt)

"
%enddef
%feature("docstring" ,tmr_chk_str) timer_check_continuation;


%feature("autodoc", "random_set_random_seed(seed)") random_set_random_seed;
%define rng_seed_str
"
    Seeds the random number generator used internally by the 
    hs_alkane library. Use seed=0 to generate seed based on 
    the current system time.

    Parameters
    ----------

    seed       : integer random number seed

"
%enddef
%feature("docstring", rng_seed_str) random_set_random_seed;

%feature("autodoc", "random_uniform_random()") random_uniform_random;
%define rng_uni_str
"
    Uses the random number generator internal to the hs_alkane
    library to return a random number on the interval [0,1).

    Returns
    ----------

    xi       : double precision random number

"
%enddef
%feature("docstring", rng_uni_str)random_uniform_random;

%feature("autodoc", "quat_axis_angle_to_quat(rotation_axis, angle)") quat_axis_angle_to_quat;
%define qt_aaq_str
"
    Returns a quaternion which represents rotation of a vector 
    around an axis by some angle.

    Parameters
    ----------
    
    rotation_axis : numpy array, unit vector about which to rotate
    angle         : rotation angle (radians)


    Returns
    ----------

    quaternion    : quaternion representing the rotation

"
%enddef
%feature("docstring", qt_aaq_str) quat_axis_angle_to_quat;

%feature("autodoc", "quat_get_minimum_arc(unit_vector1, unit_vector2)") quat_get_minimum_arc;
%define qt_mrc_str
"
    Returns a quaternion which represents rotation of unit_vector1
    onto unit_vector2 via the minimum arc

    Parameters
    ----------
    
    unit_vector1 : numpy array, unit vector before rotation
    unit_vector2 : numpy array, unit vector after rotation


    Returns
    ----------

    quaternion    : quaternion representing the rotation

"
%enddef
%feature("docstring", qt_mrc_str) quat_get_minimum_arc;

%feature("autodoc", "quat_product(quaternion1, quaternion2, normalised)") quat_product;
%define qt_prd_str
"
    Returns a quaternion which represents the rotation equivalent
    to the rotation represented by quaternion1 followed by the
    rotation represented by quaternion2. 

    Parameters
    ----------
    
    quaternion1 : numpy array, quaternion representing first rotation
    quaternion2 : numpy array, quaternion representing second rotation
    normalised  : integer flag (1/0) indicating if output quaternion
                  should be normalised


    Returns
    ----------

    quaternion    : quaternion representing combined rotation

"
%enddef
%feature("docstring", qt_prd_str) quat_product;

%feature("autodoc", "quat_inverse(quaternion)") quat_inverse;
%define qt_inv_str
"
    Returns a quaternion which represents the inverse rotation
    of a quaternion the rotation which reverses that input.

    Parameters
    ----------
    
    quaternion : numpy array, quaternion representing a rotation

    Returns
    ----------

    inverse    : numpy array, quaternion representing inverse rotation

"
%enddef
%feature("docstring", qt_inv_str) quat_inverse;

%feature("autodoc", "quat_conjugate_q_with_v(quaternion,vector)") quat_conjugate_q_with_v;
%define qt_cnj_str
"
    Conjugates a quaternion with a vector, i.e. returns a
    vector resulting from applying the rotation represented
    by the quaternion to the vector.

    Parameters
    ----------
    
    quaternion : numpy array, quaternion representing a rotation
    vector     : numpy array, vector to be rotated

    Returns
    ----------

    out_vector : numpy array, result of the conjugation

"
%enddef
%feature("docstring", qt_cnj_str) quat_conjugate_q_with_v;

%feature("autodoc", "box_set_num_boxes(num_replicas)") box_set_num_boxes;
%define bx_snb_str
"
    Sets the number of simulation boxes/replicas to be created
    and simulated. The default number of boxes is 1 and in most
    cases this should not need to be changed. 
    
    Multiple boxes may be required if implementing advanced
    sampling schemes involving multiple replicas. If so the
    number of boxes must be set before the simulation box and
    associated data structures are initialised. It cannot be
    changed during the simulation without generating
    udefined behaviour.

    Parameters
    ----------
    
    num_boxes  : number of boxes/replicas to use

"
%enddef
%feature("docstring", bx_snb_str) box_set_num_boxes;

%feature("autodoc", "box_get_num_boxes()") box_get_num_boxes;
%define bx_gnb_str
"
    Queries that the number of simulation boxes/replicas that
    are in use currently. This should not change during a simulation.

    Returns
    ----------
    
    num_boxes  : integer, number of boxes in use

"
%enddef
%feature("docstring", bx_gnb_str) box_get_num_boxes;

%feature("autodoc", "box_initialise()") box_initialise;
%define bx_init_str
"
    Allocates memory for the simulation box(es).
    
    If an initial set of cell vectors have been specified via a
    call to io_read_input() then these will be used to populate
    these data structures. Alternatively, io_real_xmol() or
    box_set_cell() can be called after box_initialise() can to
    specify the initial cell vectors for simulation.
    
    Note that cell vectors are not required if periodic boundary
    conditions are not in use - see box_set_pbc. The default is 
    that periodic boundary conditions *will* be used.

"
%enddef
%feature("docstring", bx_init_str) box_initialise;

%feature("autodoc", "box_destroy()") box_destroy;
%define bx_dsty_str
"
    Releases memory used for the simulation box(es).
    
    This should only be required if needing to reinitialise
    the simulation box(es) for a new simulation.

"
%enddef
%feature("docstring", bx_dsty_str) box_destroy;

%feature("autodoc", "box_set_cell(ibox,cell_matrix)") box_set_cell;
%define bx_scl_str
"
    Sets the cell matrix containing the three vectors
    which define the 3D simulation cell of box/replica ibox.
    
    The hs_alkane library adopts the convention that the first
    cell vector always lies along the x-axis, and the second cell 
    vector lies in the x-y plane. 
    
    Parameters
    ----------
    
    ibox         : Index (1-based) of the box/replica to set
                   cell vectors for.
    
    cell_matrix  : Numpy array compatible with that returned
                   by the get_cell() method of an ASE atoms 
                   object.

"
%enddef
%feature("docstring", bx_scl_str) box_set_cell;

%feature("autodoc", "box_get_cell(ibox)") box_set_cell;
%define bx_gcl_str
"
    Gets the cell matrix containing the three vectors
    which define the 3D simulation cell of box/replica ibox.
    
    The hs_alkane library adopts the convention that the first
    cell vector always lies along the x-axis, and the second cell 
    vector lies in the x-y plane. 
    
    Parameters
    ----------
    
    ibox         : Index (1-based) of the box/replica to get
                   cell vectors from.
                   
    Returns
    -------
    
    cell_matrix  : Numpy array compatible with that returned
                   by the get_cell() method of an ASE atoms 
                   object.    

"
%enddef
%feature("docstring", bx_gcl_str) box_get_cell;

%feature("autodoc", "box_cart_to_frac(ibox, vector)") box_cart_to_frac;
%define bx_ctf_str
"
    Converts a vector in Cartesian coordinates to one expressed
    in terms of fractional coordinates of box/replica ibox. 
    
    Parameters
    ----------
    
    ibox         : Index (1-based) of the box/replica to use
                   for conversion into fractional coordinates.
                   
    vector       : Numpy array, input Cartesian vector.
                   
    Returns
    -------
    
    out_vector   : Numpy array, output fractional vector. 
    

"
%enddef
%feature("docstring", bx_ctf_str) box_cart_to_frac;

%feature("autodoc", "box_frac_to_cart(ibox, vector)") box_frac_to_cart;
%define bx_ftc_str
"
    Converts a vector in expressed as fraction coordinates of
    box/replica ibox into the same vector in Cartesian coordinates.
    
    Parameters
    ----------
    
    ibox         : Index (1-based) of the box/replica in which
                   input fractional coordinates are expressed.
                   
    vector       : Numpy array, input fraction vector.
                   
    Returns
    -------
    
    out_vector   : Numpy array, output Cartesian vector. 
    

"
%enddef
%feature("docstring", bx_ftc_str) box_frac_to_cart;

%feature("autodoc", "box_minimum_image(ibox, vector1, vector2)") box_minimum_image;
%define bx_mic_str
"
    Finds the shortest periodic image of vector2 - vector1, i.e.
    for two position vectors finds the shortest vector which connects
    the first to the second accross a periodic boundary.
    
    Parameters
    ----------
    
    ibox         : Index (1-based) of the box/replica which defines
                   the periodic boundaries of interest.
                   
    vector1      : Numpy array, first position vector.
    
    vector2      : Numpy array, second position vector.
                   
    Returns
    -------
    
    out_vector   : Numpy array, output minimum image vector. 
    

"
%enddef
%feature("docstring", bx_mic_str) box_minimum_image;

%feature("autodoc", "box_set_pbc(flag)") box_set_pbc;
%define bx_sbc_str
"
    Can be used to set an internal variable which overrides use
    of periodic boundary conditions when calculating distances
    between position vectors in the simulation. 
    
    Should be called with flag=0 at the start of a simulation
    after box_initialise() if an open system simulation is required.
    Otherwise periodic boundary conditions will be used by default.
    
    Parameters
    ----------
    
    flag         : Integer (1 or 0) indicating if periodic
                   boundaries should be applied to the calculation
                   of distances between position vectors.  
                   
"
%enddef
%feature("docstring", bx_sbc_str) box_set_pbc;

%feature("autodoc", "box_compute_volume(ibox))") box_compute_volume;
%define bx_cvl_str
"
    Calculates the volume of the simulation cell for the specified
    simulation box/replica.
    
    Parameters
    ----------
    
    ibox         : Index (1-based) of the box/replica.
    
    Returns
    -------
    
    volume       : Volume of simulation box/replica ibox. 
                   
"
%enddef
%feature("docstring", bx_cvl_str) box_compute_volume;

%feature("autodoc", "box_set_use_verlet_list(flag)") box_set_use_verlet_list;
%define bx_svl_str
"
    If link cells are not in use then this flag can be used to
    force use of Verlet neighbour lists when computing distances
    between position vectors in the simulation box(es).
    
    * CURRENTLY NON FUNCTIONAL *
    
    Parameters
    ----------
    
    flag         : Integer (1 or 0) indicating if Verlet lists
                   should be used in the case where link cells
                   are bypassed.
    
"
%enddef
%feature("docstring", bx_svl_str) box_set_use_verlet_list;

%feature("autodoc", "box_set_link_cell_length(length)") box_set_link_cell_length;
%define bx_sll_str
"
    Sets the minimum size of a link cell in any direction in any
    simulation box. Smaller values increase the number of link
    cells and hence reduce the computational cost of finding
    distances between beads. 
    
    Link cells must be sufficiently large that all relevant
    interactions between bead i and bead j are captured within
    the 27 link cells enclosing and adjacent to that containing
    bead i. For hard sphere this means that the minimum link
    cell size in any direction must be at least a hard sphere
    diameter.
    
    The default mimimum length is 1.5 units which is appropriate
    for a bead radius of 1 provided the simulation cell angles
    remain within sensible limits, i.e. 60 to 120 degrees. 
    
    Smaller values may improve performance for simulation cells
    which are guaranteed to remain near cuboidal.
    
    Note that the subsequent use of link cells assume that any
    distortion/shinkage/expansion of the simulation box(es)
    will not invalidate this parameter.

    Parameters
    ----------
    
    length       : Minimum length of link cells in any direction.
    
"
%enddef
%feature("docstring", bx_sll_str) box_set_link_cell_length;

%feature("autodoc", "box_construct_link_cells(ibox)") box_construct_link_cells;
%define bx_clc_str
"
    Constructs and populates the link cell data structure for
    box/replica ibox. Should be called after populating 
    the initial cell_matrix of box/replica ibox.
    
    Subsequent use of link cells assumes that any
    distortion/shinkage/expansion of the simulation box(es)
    will not invalidate the ability of the link cell
    data structure to capture all relevant interactions
    between pairs of beads.   

    Parameters
    ----------
    
    ibox         : Index (1-based) of the box/replica.
    
"
%enddef
%feature("docstring", bx_clc_str) box_construct_link_cells;

%feature("autodoc", "box_destroy_link_cells(ibox)") box_destroy_link_cells;
%define bx_dlc_str
"
    Releases memory used for the link cell data structure.
    
    This should only be required if needing to reinitialise
    the link cell structure for a new simulation. 

    Parameters
    ----------
    
    ibox         : Index (1-based) of the box/replica.
    
"
%enddef
%feature("docstring", bx_dlc_str) box_destroy_link_cells;

%feature("autodoc", "box_set_bypass_link_cells(flag)") box_set_bypass_link_cells;
%define bx_blc_str
"
    Sets a flag whic instructs the hs_alkane library to bypass
    use of link cells and calculate all pairwise interactions
    directly using the minimum image convention.
    
    This can be useful if expecting the simulation cell/size
    to change significantly during simulation, or is the system
    is too small to make efficient use of link cells.
    
    The default is that links cells *will* be used.

    Parameters
    ----------
    
    flag         : Integer (1 or 0) indicating if link cells
                   should be bypassed.
    
"
%enddef
%feature("docstring", bx_blc_str) box_set_bypass_link_cells;


/* This will be parsed to generate the wrapper */
%include "timer.h"
%include "random.h"
%include "quaternion.h"
%include "box.h"
%include "alkane.h"
