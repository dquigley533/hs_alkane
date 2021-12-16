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
#include "vis_module.h"
#include "io.h"
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
%apply int *OUTPUT {int *dum_ktrial};
%apply int *OUTPUT {int *dum_max_regrow};


/* doubles */
%apply double *OUTPUT {double *rbfactor};
%apply double *OUTPUT {double *boltz_out};
%apply double *OUTPUT {double *angle};
%apply double *OUTPUT {double *dum_dr};
%apply double *OUTPUT {double *dum_dt};
%apply double *OUTPUT {double *dum_dh};
%apply double *OUTPUT {double *dum_axis};
%apply double *OUTPUT {double *dum_dv};
%apply double *OUTPUT {double *dumlength};
%apply double *OUTPUT {double *dumangle};


/* Matrices of cell vectors */
%apply(double IN_ARRAY2[ANY][ANY]) {(double cell_matrix[3][3])};
%apply(double IN_ARRAY2[ANY][ANY]) {(double delta_H[3][3])};
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

%feature("autodoc", "box_destroy_link_cells") box_destroy_link_cells;
%define bx_dlc_str
"
    Releases memory used for the link cell data structure.
    
    This should only be required if needing to reinitialise
    the link cell structure for a new simulation. 
    
"
%enddef
%feature("docstring", bx_dlc_str) box_destroy_link_cells;

%feature("autodoc", "box_set_bypass_link_cells(flag)") box_set_bypass_link_cells;
%define bx_blc_str
"
    Sets a flag which instructs the hs_alkane library to bypass
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

%feature("autodoc", "alkane_get_nchains();") alkane_get_nchains;
%define alk_gtnc_str
"
    Returns the number of chains currently in use within the simulation.


    Returns
    -------
    
    nchains      : Integer number of chains. 

"
%enddef
%feature("docstring", alk_gtnc_str) alkane_get_nchains;

%feature("autodoc", "alkane_set_nchains(nchains);") alkane_set_nchains;
%define alk_stnc_str
"

    Sets the number of chains to use within the simulation. Subsequent
    calls to alkane_init will allocate data structures to hold this 
    number of chains. 

    Can be modified during the simulation. This is useful for when the
    number of chains is changing, such as when growing new random 
    chain configurations to avoid checking for overlaps with chains
    that have not yet been created.

    Parameters
    ----------

    nchains      : Integer number of chains to use during the simulation. 


"
%enddef
%feature("docstring", alk_stnc_str) alkane_set_nchains;

%feature("autodoc", "alkane_get_nbeads(int_*nbeads);") alkane_get_nbeads;
%define alk_gtnb_str
"

    Returns the number of beads per chain currently in use within the 
    simulation.

    Returns
    -------
    
    nbeads      : Integer number of beads per chain 

"
%enddef
%feature("docstring", alk_gtnb_str) alkane_get_nbeads;

%feature("autodoc", "alkane_set_nbeads(nbeads);") alkane_set_nbeads;
%define alk_stnb_str
"

    Sets the number of beans per chain to use within the simulation. 
    Subsequent calls to alkane_init will allocate data structures to hold
    this many beads per chain. 

    Parameters
    ----------

    nbeads      : Bumber of beads per chain to use during the simulation. 

"
%enddef
%feature("docstring", alk_stnb_str) alkane_set_nbeads;

%feature("autodoc", "alkane_initialise();") alkane_initialise;
%define alk_init_str
"

    Creates data structures to hold information on the current configuration
    of beads/chains in all simualtion boxes/replicas. Must be called before
    calling routines which implement trial Monte Carlo moves on the configuration.

    Requires that the simulation box has been initialised via box_initialise().

    This routine can be invoked multiple times within a session. Each call
    will result in destruction of the existing data structures and their
    contents before creating new ones. This can be used to implement
    multiple simulations within one Python session.


"
%enddef
%feature("docstring", alk_init_str) alkane_initialise;

%feature("autodoc", "alkane_construct_linked_lists(int_ibox);") alkane_construct_linked_lists;
%define alk_cnstrctll_str
"
    Creates and populates the linked list data structures used to 
    accelerate computation of neighbour distances in box ibox.
    
    Can be invoked manually to ensure current linked lists are 
    valid, but is normally called only by other routine which
    change the size/shape of the simulation cell.

    See alkane_update_linked_lists for single bead updates of the 
    linked list data structures. 

    Link cell data structures in the box module must be up-to-date
    when calling this rouine. See box_construct_link_cells().

    Has no effect if use of link cells has been bypassed see
    box_set_bypass_link_cells().

    Parameters
    ----------

    ibox        : Box/replica for which to construct linked lists

"
%enddef
%feature("docstring", alk_cnstrctll_str) alkane_construct_linked_lists;

%feature("autodoc", "alkane_destroy();") alkane_destroy;
%define alk_dstry_str
"
    Destroys all data structures related to beads/chains in the
    simulation cell. Acts on all boxes/replicas.

    Use alkane_initialise() to both delete these data structures
    and create new ones.

"
%enddef
%feature("docstring", alk_dstry_str) alkane_destroy;

%feature("autodoc", "alkane_grow_chain(ichain,ibox,new_conf);") alkane_grow_chain;
%define alk_grwc_str
"
    Grows a new chain configuration for chain ichain in box 
    ibox from a randomly selected bead onwards. Uses Rosenbluth
    sampling. Primarily used for configurational bias Monte Carlo.

    Invoke with new_conf=0 to compute the Rosenbluth factor for 
    the existing chain configuration, and new_conf=0 to compute
    the Rosenbluth factor for a new trial configuration

    If the chain has not previously been assigned coordinates then
    a new chain is grown from the first bead if new_conf=1. This 
    can be useful in creating random initial configurations. 

    If ifail is non-zero then the algorithm has failed to generate
    a trial chain segment which has no overlaps with other beads. 
    The number of trial attempts (ktrial) can be changed via
    alkane_set_ktrial().

    The maximum number of beads which are (re)grown by the trial
    move is controlled via alkane_set_max_regrow().


    Parameters
    ----------

    ichain     :  Integer (1-based) chain to grow/regrow
    ibox       :  Integer (1-based) box/replica to use
    new_conf   :  1 to grow new chain, 0 otherwise

    Returns
    -------

    rb_factor  : Rosenbluth factor for existing/trial chain
    ifail      : non-zero integer if chain growth failed


"
%enddef
%feature("docstring", alk_grwc_str) alkane_grow_chain;

%feature("autodoc", "alkane_get_chain(ichain,ibox);") alkane_get_chain;
%define alk_gtch_str
"
    Returns a Numpy array containing the bead coordinates of chain 
    ichain in box/replica ibox.

    NOTE : Note the returned Numpy array is a 'view' of the internal
    representation of coordinates. Updating the elements of the Numpy
    array will update the coordinates inside the model. This avoids
    making copies. i.e. 

    mychain = mdl.alkane_get_chain(ichain, ibox)
    mychain[0][2] = 3.4

    will modify the z coordinate of bead zero in the chain within
    the model.

    Note also however that assigning the name of the Numpy array  
    to a new array will *not* update the model, so

    mychain = mdl.alkane_get_chain(ichain, ibox)
    mychain = np.zeros(3,Nbeads)

    will not set all coordinate for the chain to zero. Instead
    it will replace the Numpy array holding the 'view' of the 
    chain coordinates with another empty array not linked to 
    the internal model data.

    Parameters
    ----------

    ichain     :  Integer (1-based) chain index
    ibox       :  Integer (1-based) box/replica index

    Returns
    -------
    
    chain_out  :  Numpy array with 'view' of model chain coordinates

"
%enddef
%feature("docstring", alk_gtch_str) alkane_get_chain;

%feature("autodoc", "alkane_translate_chain(ichain,ibox);") alkane_translate_chain;
%define alk_trch_str
"
    Makes a trial translation move by translating an entire chain in
    a random direction by a distance between zero and dr_max.

    Returns the probability of accepting the move, which will be 1
    if the move generates no overlaps with other beads, zero otherwise.
    
    See also alkane_set_dr_max() and alkane_get_dr_max().

    Parameters
    ----------

    ichain     :  Integer (1-based) chain index to translate
    ibox       :  Integer (1-based) box/replica index to use


    Returns
    -------
    
    boltz      :  Probability of accepting the trial move

"
%enddef
%feature("docstring", alk_trch_str) alkane_translate_chain;

%feature("autodoc", "alkane_rotate_chain(ichain,ibox,bond);") alkane_rotate_chain;
%define alk_rtch_str
"
    Makes a trial rotation move by rotating an entire chain about
    a random axis by an angle between zero and dt_mx.

    If the input parameter bond is equal to 1 then the chain will
    be rotated around an axis defined by the vector connecting
    bead 0 to bead 1 on that chain. For dense packed systems of linear 
    chains replacing some (large) fraction of standard rotation 
    moves with moves of this type can improve sampling efficiency. 
    In this case the random rotation is between zero and axis_max.

    Returns the probability of accepting the move, which will be 1
    if the move generates no overlaps with other beads, zero otherwise.
    The quaternion representing the trial rotation is also returned.

    See also alkane_set_dr_max() and alkane_get_dr_max() or
             alkane_set_axis_max() and alkane_get_axis_max()
  
    Parameters
    ----------

    ichain     :  Integer (1-based) chain index to rotate
    ibox       :  Integer (1-based) box/replica index to use
    bond       :  Integer (1 or 0) rotate around first bond

    Returns
    -------
    
    boltz      :  Probability of accepting the trial move
    quat       :  Quaternion representing the rotation 

"
%enddef
%feature("docstring", alk_rtch_str) alkane_rotate_chain;

%feature("autodoc", "alkane_bond_rotate(ichain,ibox,allow_flip);") alkane_bond_rotate;
%define alk_dhmv_str
"
    Makes a trial move in which a randomly selected dihedral angle
    on chain ichain in box/replica ibox is rotated by an angle
    between zero and dh_max. Additionally, if allow_flip=1
    there is a probabilty of 50% that an additional 
    2Pi/3 radians is added such that the chains flips between
    'gauche' and 'anti' conformations about this bond.

    See also alkane_set_dh_max() and alkane_get_dh_max()

    Parameters
    ----------

    ichain     :  Integer (1-based) chain index to use
    ibox       :  Integer (1-based) box/replica index to use
    allow_flip :  Integer (1 or 0) include flips between conformations

    Returns
    -------

    ia         : First bead involved in the dihedral angle
    angle      : Angle of rotation about the bond
    boltz      : Probability of accepting the trial move

"
%enddef
%feature("docstring", alk_dhmv_str) alkane_bond_rotate;

%feature("autodoc", "alkane_check_chain_overlap(ibox);") alkane_check_chain_overlap;
%define alk_chkov_str
"
    Checks if simulation box/replica ibox contains any overlaps
    between beads. Both inter and intra-chain overlaps are
    detected. As specified by the choice of model, interactions
    between beads on the same chain seperated by 1, 2 and (in 
    some cases) 3 bonds are ignored.
  
    Useful as a sanity check. There should never be any overlaps
    if all list structures are up-to-date and no moves which 
    generate overlaps have been accepted.

    Parameters
    ----------

    ibox       :  Integer (1-based) box/replica index to use
    
    Returns
    -------

    overlap    :  Integer (1/0) indicating if overlaps are present

"
%enddef
%feature("docstring", alk_chkov_str) alkane_check_chain_overlap;

%feature("autodoc", "alkane_check_chain_geometry(int_ichain,_int_ibox,_int_*violated);") alkane_check_chain_geometry;
%define alk_chkgm_str
"
    Checks that the internal geometry of chain ichain in box/replica
    ibox is consistent with the model bond length and bond angle  
    constraints.

    Parameters
    ----------

    ichain     :  Integer, chain to check (1-based)
    ibox       :  Integer (1-based) box/replica index to use

    Returns
    -------

    violated   :  Integer (1/0) indicatign if contraints are violated

"
%enddef
%feature("docstring", alk_chkgm_str) alkane_check_chain_geometry;

%feature("autodoc", "alkane_update_linked_lists(ibead,ichain,ibox,old_pos[3],new_pos[3]);") alkane_update_linked_lists;
%define alk_updtll_str
"
    Updates the linked list data structure to account for bead ibead
    on chain ichain in box ibox having moved from position old_pos
    to new_pos. 

    In principle this should be used after every accepted trial move
    which changes a bead position. However for mostly static solids, 
    the set of link cells which contains all neighbours of each bead
    is unlikely to change and link list updates can (with caution) be
    omitted.

 
    Parameters
    ----------

    ibead      :  Index (1-based) of bead on ichain to update
    ichain     :  Integer, chain on which bead resides (1-based)
    ibox       :  Integer (1-based) box/replica index to use
    old_pos    :  Numpy array - old x,y,z coords of bead
    new_pos    :  Numpy array - new x,y,z coords of bead

"
%enddef
%feature("docstring", alk_updtll_str) alkane_update_linked_lists;

%feature("autodoc", "alkane_construct_neighbour_list(ibox);") alkane_construct_neighbour_list;
%define alk_cnsctnl_str
"
    Creates a Verlet neighbour for all beads in simulation box/replica ibox.
    A neighbour list can be used instead of link cells to accelerate
    checks for overlaps between builds. See

    box_set_bypass_link_cells()
    box_set_use_verlet_list()

    All neighbours of a bead within a distance of two bead diameters
    are added into the neighbour list.

    There is no automatic update of the neighbour list to account for 
    movement of beads. This function should be called whenever 
    neighbours are likely to have moved sufficiently that the previous
    neighbour list has been invalidatd.

 
    Parameters
    ----------

    ibox       :  Integer (1-based) box/replica index to construct for

"
%enddef
%feature("docstring", alk_cnsctnl_str) alkane_construct_neighbour_list;

%feature("autodoc", "alkane_get_dr_max();") alkane_get_dr_max;
%define alk_gtdrm_str
"
    Queries the maximum distance a chain will be moved during a 
    trial Monte Carlo translation.

    Returns
    -------
    
    dr_max     :  Maximum displacement

"
%enddef
%feature("docstring", alk_gtdrm_str) alkane_get_dr_max;

%feature("autodoc", "alkane_set_dr_max(dr_max);") alkane_set_dr_max;
%define alk_stdrm_str
"
    Sets the maximum distance a chain will be moved during a 
    trial Monte Carlo translation.

    Parameters
    ----------
    
    dr_max     :  Maximum displacement

"
%enddef
%feature("docstring", alk_stdrm_str) alkane_set_dr_max;

%feature("autodoc", "alkane_get_dt_max();") alkane_get_dt_max;
%define alk_gtdtm_str
"
    Queries the maximum angle a chain will be rotated around
    a random axis during a Monte Carlo trial rotation.

    Returns
    -------
    
    dt_max     :  Maximum rotation angle


"
%enddef
%feature("docstring", alk_gtdrm_str) alkane_get_dt_max;

%feature("autodoc", "alkane_set_dt_max(dt_max);") alkane_set_dt_max;
%define alk_stdtm_str
"
    Sets the maximum angle a chain will be rotated around
    a random axis during a Monte Carlo trial rotation.

    Parameters
    ----------
    
    dt_max     :  Maximum rotation angle



"
%enddef
%feature("docstring", alk_stdtm_str) alkane_set_dt_max;

%feature("autodoc", "alkane_get_axis_max();") alkane_get_axis_max;
%define alk_gtaxm_str
"
    Queries the maximum angle a chain will be rotated around
    an axis through its first bond vector during a Monte
    Carlo trial rotation.

    Returns
    -------
    
    axis_max     :  Maximum rotation angle around bond axis

"
%enddef
%feature("docstring", alk_gtaxm_str) alkane_get_axis_max;

%feature("autodoc", "alkane_set_axis_max(axis_max);") alkane_set_axis_max;
%define alk_staxm_str
"
    Sets the maximum angle a chain will be rotated around
    an axis through its first bond vector during a Monte
    Carlo trial rotation.

    Parameters
    ----------
    
    axis_max     :  Maximum rotation angle around bond axis

"
%enddef
%feature("docstring", alk_staxm_str) alkane_set_axis_max;

%feature("autodoc", "alkane_get_dh_max();") alkane_get_dh_max;
%define alk_gtdhm_str
"
    Queries the maximum change in dihedral angle made during
    a trial Monte Carlo rotation about a random bond.

    Returns
    -------
    
    dh_max       :  Maximum rotation angle around bond 

"
%enddef
%feature("docstring", alk_gtdhm_str) alkane_get_dh_max;

%feature("autodoc", "alkane_set_dh_max(dh_max);") alkane_set_dh_max;
%define alk_stdhm_str
"
    Sets the maximum change in dihedral angle made during
    a trial Monte Carlo rotation about a random bond.

    Returns
    -------
    
    dh_max       :  Maximum rotation angle around bond 

"
%enddef
%feature("docstring", alk_stdhm_str) alkane_set_dh_max;

%feature("autodoc", "alkane_get_dv_max();") alkane_get_dv_max;
%define alk_gtdvm_str
"
    Queries the parameter used to determine the maximum 
    change in the simulation cell during a trial move as
    implemented in alkane_box_resize().

    For anisotopic box moves this parameter controls the
    maximum change in any component of a randomly selected
    cell vector. For isotropic (pure expansion/contraction)
    volume moves the parameter controls the maximum 
    change in volume.

    See also box_get_isotropic().

    Returns
    -------
    
    dv_max       :  Maximum box change parameter
"
%enddef
%feature("docstring", alk_gtdvm_str) alkane_get_dv_max;

%feature("autodoc", "alkane_set_dv_max(dv_max);") alkane_set_dv_max;
%define alk_stdvm_str
"
    Sets the parameter used to determine the maximum 
    change in the simulation cell during a trial move as
    implemented in alkane_box_resize().

    For anisotopic box moves this parameter controls the
    maximum change in any component of a randomly selected
    cell vector. For isotropic (pure expansion/contraction)
    volume moves the parameter controls the maximum 
    change in volume.

    See also box_set_isotropic().

    Parameters
    ----------
    
    dv_max       :  Maximum box change parameter

"
%enddef
%feature("docstring", alk_stdvm_str) alkane_set_dv_max;

%feature("autodoc", "alkane_get_ktrial();") alkane_get_ktrial;
%define alk_gtktr_str
"
    Queries the parameter which controls the number of trial
    chain segments used at each bead when (re)growing chains
    in alkane_grow_chain(). 

    Returns
    -------

    ktrial      :   Number of trials  

"
%enddef
%feature("docstring", alk_gtktr_str) alkane_get_ktrial;

%feature("autodoc", "alkane_set_ktrial(int_ktrial);") alkane_set_ktrial;
%define alk_stktr_str
"
    Sets the parameter which controls the number of trial
    chain segments used at each bead when (re)growing chains
    in alkane_grow_chain(). 

    Parameters
    ----------

    ktrial      :   Number of trials  

"
%enddef
%feature("docstring", alk_stktr_str) alkane_set_ktrial;

%feature("autodoc", "alkane_get_max_regrow();") alkane_get_max_regrow;
%define alk_gtmxrg_str
"
    Queries the maximum number of beads/segments to (re)grow 
    during a configurational bias Monte Carlo move as 
    implemented in alkane_grow_chain().

    Returns
    -------

    max_regrow  : Maximum number of beads/segments to regrow

"
%enddef
%feature("docstring", alk_gtmxrg_str) alkane_get_max_regrow;

%feature("autodoc", "alkane_set_max_regrow(int_dum_max_regrow);") alkane_set_max_regrow;
%define alk_stmxrg_str
"
    Sets the maximum number of beads/segments to (re)grow 
    during a configurational bias Monte Carlo move as 
    implemented in alkane_grow_chain().

    Parameters
    ----------

    max_regrow  : Maximum number of beads/segments to regrow

"
%enddef
%feature("docstring", alk_stmxrg_str) alkane_set_max_regrow;

%feature("autodoc", "alkane_change_box(ibox,delta_H[3][3]);") alkane_change_box;
%define alk_chbx_str
"
    Changes the matrix of simulation cell vectors in box ibox
    by a matrix delta_H (3x3 Numpy array). Fractional positions
    of the first bead in each chain are preserved.
 
    This may lead to chain overlaps which are not tested. Note that
    this function does not implement a trial move to be accepted
    or rejected. It is provided for the purposes of implementing
    higher-level sampling algorithms in the calling routine.

    Parameters
    ----------

    ibox        : Simulation box/replica to change
    delta_H     : 3x3 Numpy array - change in matrix of cell vectors

"
%enddef
%feature("docstring", alk_chbx_str) alkane_change_box;

%feature("autodoc", "alkane_box_scale(int_ibox,_double_scaleA,_double_scaleB,_double_scaleC);") alkane_box_scale;
%define alk_bxscl_str
"
    Changes the matrix of simulation cell vectors in box ibox
    by scaling the three cell vectors by the factors specified.
    Fractional positions of the first bead in each chain are preserved.

    This may lead to chain overlaps which are not tested. Note that
    this function does not implement a trial move to be accepted
    or rejected. It is provided for the purposes of implementing
    higher-level sampling algorithms in the calling routine.

    Parameters
    ----------

    ibox        : Simulation box/replica to change
    scaleA      : Scaling factor for first  cell vector
    scaleB      : Scaling factor for second cell vector
    scaleC      : Scaling factor for third  cell vector 

"
%enddef
%feature("docstring", alk_bxscl_str) alkane_box_scale;

%feature("autodoc", "alkane_box_resize(pressure,ibox,reset);") alkane_box_resize;
%define alk_bxrsz_str
"
    Implements a trial Monte Carlo move in which the size/shape
    of simulation box is changed by a random amount controled
    by the parameter dv_max.

    For isotropic moves a random change involume between 
    -dv_max and +dv_max is proposed. For anisotropic moves a
    randomly selected component of cell vector is changed by
    an amount between -dv_max and +dv_max.

    Fractional positions of the first bead in each chain are
    preserved in each case.

    See also alkane_set_dv_max(), box_set_isotropic().
   
    If the parameter reset=1 then the most recently proposed
    move of this kind is reversed, implementing rejection.


    Parameters
    ----------

    Pressure     : Simulation pressure in reduced units
    ibox         : Simulation box/replica to modify
    reset        : Integer (0/1) reverse previous box resize move

    Returns
    -------

    boltz        : Probability of accepting move

"
%enddef
%feature("docstring", alk_bxrsz_str) alkane_box_resize;

%feature("autodoc", "write_psf(nbeads,nchains);") write_psf;
%define vis_wrtpsf_str
"
    Write a Protein Structure (topology) files to chain.psf. This file
    contains bonding information and can be used in the VMD molecular
    simulation software to definate a structure into which coordinates
    are then read from dcd file(s). 
 
    If the simulation consists of multiple boxes/replicas then one file
    is created for each, i.e. chain.psf.01 chain.psf.02 etc.

    See also write_dcd_header(), write_dcd_snapshot().

    Parameter
    ---------

    nbeads    : Number of beads per chain in the simulation(s)
    nchains   : Number of chains in the simulation(s)

"
%enddef
%feature("docstring", vis_wrtpsf_str) write_psf;

%feature("autodoc", "write_dcd_header(nbeads,nchains);") write_dcd_header;
%define vis_wrtdcdhd_str
"
    Writes the header of a Charmm-style dcd file to which snapshots
    of bead coordinates can subsequently be dumped via write_dcd_snapshot.
    The header is written to chain.dcd, which is overwritten if it 
    already exists.

    If the simulation consists of multiple boxes/replicas then one file
    is created for each, i.e. chain.dcd.01 chain.dcd.02 etc.

    See also write_psf(), write_dcd_snapshot().

    Parameter
    ---------

    nbeads    : Number of beads per chain in the simulation(s)
    nchains   : Number of chains in the simulation(s)

"
%enddef
%feature("docstring", vis_wrtdcdhd_str) write_dcd_header;

%feature("autodoc", "write_dcd_snapshot();") write_dcd_snapshot;
%define vis_wrtdcdsn_str
"
    Appends a snapshot of the current model configuration to the chain.dcd
    file previously created by a call to write_dcd_header().

    If the simulation consists of multiple boxes/replicas then one file
    is used for each, i.e. chain.dcd.01 chain.dcd.02 etc.

    The dcd file(s) can be read by the VMD visualisation software 
    by loading them into a configuration previous read from chain.psf.

    See also write_psf(), write_dcd_header().
"
%enddef
%feature("docstring", vis_wrtdcdsn_str) write_dcd_snapshot;

%feature("autodoc", "io_read_xmol(filename);") io_read_xmol;
%define io_rdxmol_str
"
    Reads an initial configuration for the simulation(s) from
    filename in the current working directory. This should be 
    structured as a standard xyz file, with the 9 components
    of the matrix of cell vectors listed on the second line.

    The nbeads beads on the first chain should appear sequentially
    as entries 1-nbeads, followed by the second chain, and so on.

    The data structures inside the alkane module must have been
    created via alkane_initialise() before calling this function.

    If the simulation consists of multiple boxes/replicas then one file
    is read for each, i.e. chain.xmol.0001 chain.xmol.0002 etc.

    Parameter
    ---------

    filename   : Filename to read (optional). Uses chain.xmol otherwise.
 


"
%enddef
%feature("docstring", io_rdxmol_str) io_read_xmol;

%feature("autodoc", "io_write_xmol(filename);") io_write_xmol;
%define io_wtxmol_str
"
    Write current configuration for the simulation(s) to
    filename in the current working directory. This will be 
    structured as a standard xyz file, with the 9 components
    of the matrix of cell vectors listed on the second line.

    The data structures inside the alkane module must have been
    created via alkane_initialise() and populated with data
    efore calling this function.

    If the simulation consists of multiple boxes/replicas then one file
    is written for each, i.e. final.xmol.0001 final.xmol.0002 etc.

    Parameter
    ---------

    filename   : Filename to read (optional). Uses final.xmol otherwise.

"
%enddef
%feature("docstring", io_wtxmol_str) io_write_xmol;


%feature("autodoc", "alkane_set_bondlength(L);") alkane_set_bondlength;
%define alk_stbl_str
"

    Sets the distance between adjacent beads in a chain in a simulation.
    The bond length is expected to be in terms of the diameter of the sphere, such that
    two beads which are tangentially touching will have a bond length of 1.0. 
    

    Parameters
    ----------

    bondlength      : Floating point number to use


"
%enddef
%feature("docstring", alk_stbl_str) alkane_set_bondlength;


%feature("autodoc", "alkane_get_bondlength();") alkane_get_bondlength;
%define alk_gtbl_str
"
    Returns the bond length currently in use within the simulation.


    Returns
    -------
    
    L      : Distance between the centres of two adjacent beads in a chain. 

"
%enddef
%feature("docstring", alk_gtbl_str) alkane_get_bondlength;


%feature("autodoc", "alkane_set_bondangle(bondangle);") alkane_set_bondangle;
%define alk_stba_str
"

    Sets the angle formed by three adjacent beads in a chain in a simulation, in degrees
. 
    

    Parameters
    ----------

    bondangle      : Angle to set between beads in degrees


"
%enddef
%feature("docstring", alk_stba_str) alkane_set_bondangle;


%feature("autodoc", "alkane_get_bondangle();") alkane_get_bondangle;
%define alk_gtba_str
"
    Returns the angle formed by three adjacent beads in a chain in a simulation, in degrees

    Returns
    -------
    
    bondangle      : Angle between the centres of three consecutive beads in a chain. 

"
%enddef
%feature("docstring", alk_gtba_str) alkane_get_bondangle;



/* This will be parsed to generate the wrapper */
%include "timer.h"
%include "random.h"
%include "quaternion.h"
%include "box.h"
%include "alkane.h"
%include "vis_module.h"

/* Explicit here rather than included due to optional argument */
void io_read_xmol(char* filename="chain.xmol");
void io_write_xmol(char* filename="final.xmol");
