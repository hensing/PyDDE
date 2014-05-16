#!/usr/bin/env python
"setup for PyDDE"
import sys

from warnings import warn
from setuptools import setup, find_packages, Extension

try:
    from numpy import get_include
    INCLUDE_DIRS = [get_include()]
except ImportError:
    warn("numpy not found!")
    INCLUDE_DIRS = []


REQUIREMENTS = [
    'numpy',
]

# mac osx 10.9 clang fix
if sys.platform == 'darwin':
    EXTRA_COMPILE_ARGS = [
        '-Wno-error=unused-command-line-argument-hard-error-in-future']
else:
    EXTRA_COMPILE_ARGS = []


EXT_MODULES = [
    Extension(
        "PyDDE.ddesolve",
        sources=['PyDDE/src/ddeq.c', 'PyDDE/src/ddesolve95.c',
                 'PyDDE/src/wrapper.c'],
        extra_compile_args = EXTRA_COMPILE_ARGS),
]


setup(
    name='PyDDE',
    version='0.2.2',
    description="PyDDE is a solver for delay differential equations \
                 written in Python and C.",
    long_description="""
    PyDDE is an open source numerical solver for systems of delay differential
    equations (DDEs), implemented as a Python package and written in both
    Python and C.  It is built around the numerical routines of the R package
    ddesolve, which is itself based on Simon Wood's Solv95, a DDE solver for
    Microsoft Windows systems written in C.
    """,
    classifiers=[
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
    keywords="delay differential equation solver dde switches solv95 ddesolve",
    maintainer="Henning Dickten",
    maintainer_email="hdickten@uni-bonn.de",
    url="https://www.github.com/hensing/PyDDE",
    #py_modules=['PyDDE.pydde'],
    requires=REQUIREMENTS,
    #extras_require=extras_require,
    packages=find_packages(),
    #packages=['PyDDE'],
    #scripts=scripts,
    ext_modules=EXT_MODULES,
    include_package_data=True,
    include_dirs=INCLUDE_DIRS,
    package_data={'PyDDE' : ['doc/*.pdf', 'test/*.py']},
)
