Please see the file INSTALL.txt for information on installing PyDDE.

The file LICENSE.txt for important licensing and warranty information; the 
file GPL.txt contains the actual license.  Your use of PyDDE and its source 
code is governed by these conditions; if you do not agree with them please 
uninstall PyDDE and delete the source distribution from your computer.


ABOUT PyDDE
===========

PyDDE is an open source numerical solver for systems of delay differential 
equations (DDEs), implemented as a Python package and written in both Python 
and C.  It is built around the numerical routines of the R package ddesolve, 
which is itself based on Simon Wood's Solv95, a DDE solver for Microsoft 
Windows systems written in C.


RELEASE NOTES
=============

0.2.2  -  07 December 2008

-- Fixed a bug that resulted in failure to compile on some Mac OS X systems.
   (Thanks to Josh Lippai who reported the bug and coded the fix.)
-- Various minor typographical corrections to the code and manual.

0.2.1  -  31 October 2007

-- Migrated to use the back-end from the R package 'ddesolve'.

0.1.3  -  25 April 2007

-- Fixed a memory 'leak' due to poor reference handling of Python objects in the 
   C code.  Python was keeping a reference to objects that should have been 
   deleted.
-- Fixed a segfault when updating or cleaning up the history buffer for models 
   without any history variables.
-- Fixed an apparent bug in setup.py.  The libraries field doesn't work with 
   some versions of Python.  Removing it doesn't seem to hurt, either. 

0.1.2  -  23 April 2007

-- Improved memory management by better allocating and freeing the simulated
   data and history buffer at each dde.solve() call.

0.1.1  -  13 April 2007

-- Fixed a bug whereby a data variable would be defined twice, causing 
   compilation errors on some systems (e.g. OS X with gcc 4.0.1).

0.1.0  -  14 December 2005

-- First release of PyDDE.

