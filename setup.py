##############################
# Setup script for hs_alkane #
##############################
#!/usr/bin/env python

def configuration():

    from numpy.distutils.misc_util import Configuration

    config = Configuration('hs_alkane', parent_name=None, top_path=None)

    alkane_src = ['src/constants.f90',
                  'src/timer.f90',
                  'src/random.f90',
                  'src/quaternion.F90',
                  'src/box.f90',
                  'src/alkane.f90',
                  'src/vis_module.f90',
                  'src/mc_dummy.f90',
                  'src/io.f90']
    
    config.add_library('alkane', sources=alkane_src,
#                       extra_f90_compile_args=['-g','-fbounds-check','-fbacktrace']
#                       extra_f90_compile_args=['-O3']
                       )

    config.add_extension('_alkane',
        sources = ['src/alkane.i',
                   'include/timer.h',
                   'include/random.h',
                   'include/box.h',
                   'include/quaternion.h',
                   'include/alkane.h'],
        libraries = ['alkane'],
        include_dirs = ['./include'],
        depends = ['src/alkane.i'],
    )

    config.version = "0.1.1"

    return config

if __name__ == "__main__":

    from numpy.distutils.core import setup

    setup(configuration=configuration, zip_safe=True)



