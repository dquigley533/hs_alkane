# meson.build for hs_alkane
# D. Quigley - University of Warwick

project(
	'hs_alkane',       # The name of the project
	['fortran', 'c'],  # The languages that Meson should initialize

	# Keyword arguments:
	version : '0.1.1',
	default_options : ['optimization=2', 'debug=false', 'warning_level=3'],
)

# List of options to pass to fortran compiler
# add_project_arguments(['-fbacktrace', '-fbounds-check'], language:'fortran') 

# Use swig to generate wrapper code
swig = find_program('swig')

# Need Python headers to compile the wrapper code
pydep = dependency('python3', version : '>=3.0.0')
py = import('python').find_installation(pure: false)

# We'll link using gcc so need to specify some fortran dependencies
fc = meson.get_compiler('fortran')
libgfortran_dep = fc.find_library('gfortran', required: true)
libquadmath_dep = fc.find_library('quadmath', required: true)

# I have issues on MacOS if the Apple-provided clang compiler is used. 
cc = meson.get_compiler('c')
if cc.get_id() != fc.get_id()
  message('Selected C compiler does not have same ID as Fortran compiler.')
  message('This can cause problems with linking.')
  message('Set CC and FC environment variables to compatible compilers before building')
  error()
endif


# NumPy include directory - needed in all submodules
incdir_numpy = run_command(py,
  [
    '-c',
    'import numpy; print(numpy.get_include())'
  ],
  check: false
).stdout().strip()

# List of include directories
incdir = include_directories('include')

wrapper_code = custom_target('wrapper_code',
	     input : ['src/alkane.i'],
	     output : ['alkane_wrap.c', 'alkane.py'],
	     command : [swig, '-python', '-I../include/', '-outcurrentdir', '@INPUT@'],
	     )

# Source files for fortran library
alkane_src = ['src/constants.f90',
                      'src/timer.f90',
              	      'src/random.f90',
                      'src/quaternion.F90',
                      'src/box.f90',
                      'src/alkane.f90',
                      'src/vis_module.f90',
                      'src/mc_dummy.f90',
                      'src/io.f90']
		      

# Build fortran source and wrapper code into extension module
py.extension_module(
  '_alkane',
  [alkane_src, wrapper_code],
  install: true,
  include_directories : [incdir, incdir_numpy],
  subdir: 'hs_alkane',
  dependencies: [libgfortran_dep, libquadmath_dep],

)

# Finally define and install python sources
python_sources = [
  '__init__.py',
]
py.install_sources(python_sources, subdir : 'hs_alkane')
