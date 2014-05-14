from distutils.core import setup, Extension
from distutils.sysconfig import get_config_var
setup(name="PyDDE", 
      version="0.2.2",
      description="PyDDE is a solver for delay differential equations written in Python and C.",
      long_description=""" PyDDE is an open source numerical solver for systems of delay differential equations (DDEs), implemented as a Python package and written in both Python and C.  It is built around the numerical routines of the R package ddesolve, which is itself based on Simon Wood's Solv95, a DDE solver for Microsoft Windows systems written in C.  

PyDDE can solve a wide range of ODE and DDE models with discontinuities that may have state-dependent effects but state-independent timings.  Simulation is handled by an adaptively-stepping embedded RK2(3) scheme with cubic Hermite interpolation for calculation of delay terms.  Some of the advantages of PyDDE are that it is fast, efficient and allows rapid prototyping of scriptable models in a free, platform-independent environment.  """,
      classifiers = [
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License (GPL)',
          'Operating System :: OS Independent',
          'Topic :: Scientific/Engineering :: Mathematics',
          'Topic :: Software Development :: Libraries :: Python Modules',
          ],
      platforms="Any",
      license="GPL",
      keywords="delay differential equation solver dde switches solv95 ddesolve pydde",
      maintainer="Benjamin J. Cairns",
      maintainer_email="ben.cairns@ceu.ox.ac.uk",
      url="http://users.ox.ac.uk/~clme1073/python/PyDDE/",
      py_modules=['PyDDE.pydde'],
      ext_modules=[Extension("PyDDE.ddesolve", 
                             ["PyDDE/src/ddeq.c", "PyDDE/src/ddesolve95.c", "PyDDE/src/wrapper.c"],
                             # Not sure why, but the next line doesn't work
                             # with some Python versions. Removing it is OK.
                             #libraries=["python"+get_config_var('VERSION')],
                             extra_compile_args=["-I."]
                             )
                  ],
      packages=['PyDDE'],
      package_data={'PyDDE' : ['doc/*.pdf'], 'PyDDE' : ['test/*.py']}
      )
