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
%}

/* Numpy array typemap */
%include "numpy.i"

%init %{
  import_array();
%}

/* Map array array onto arguments used in the C interface */

/* Vectors */
%apply(double IN_ARRAY1[ANY]) {(double unit_vector1[3])};
%apply(double IN_ARRAY1[ANY]) {(double unit_vector2[3])};
%apply(double IN_ARRAY1[ANY]) {(double rotation_axis[3])};
%apply(double IN_ARRAY1[ANY]) {(double vector[3])};
%apply(double ARGOUT_ARRAY1[ANY]) {(double vec_out[3])};


/* Quaternions */
%apply(double ARGOUT_ARRAY1[ANY]){(double quat[4])};
%apply(double IN_ARRAY1[ANY]){(double quaternion[4])};
%apply(double IN_ARRAY1[ANY]){(double quaternion1[4])};
%apply(double IN_ARRAY1[ANY]){(double quaternion2[4])};
%apply(double ARGOUT_ARRAY1[ANY]){(double quat_product[4])};
%apply(double ARGOUT_ARRAY1[ANY]){(double quat_inverse[4])};

/* Number of boxes */
%apply(int ARGOUT_ARRAY1[ANY]){(int num_replicas[1])};

/* Matrices of cell vectors */
%apply(double IN_ARRAY2[ANY][ANY]) {(double cell_matrix[3][3])};
%apply(double ARGOUT_ARRAY2[ANY][ANY]) {(double outmat[3][3])};


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


/* This will be parsed to generate the wrapper */
%include "timer.h"
%include "random.h"
%include "quaternion.h"
%include "box.h"
