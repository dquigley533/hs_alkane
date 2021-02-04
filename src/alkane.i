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
%}

/* Numpy array typemap */
%include "numpy.i"

%init %{
  import_array();
%}

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


%feature("autodoc", "random_set_random_seed()") random_set_random_seed;
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



/* This will be parsed to generate the wrapper */
%include "timer.h"
%include "random.h"
