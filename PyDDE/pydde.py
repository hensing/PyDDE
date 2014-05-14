'''
PyDDE -- A Python port of the C DDE solver, Solv95.

Original C code:  Simon Wood <s.wood@bath.ac.uk>
Python/C port:    Benjamin Cairns <ben.cairns@ceu.ox.ac.uk>
Last updated:     07 December 2008

This module primarily provides the class pydde.dde.

The class pydde.dde has the following methods:
    
  <instance>.initsolver   -- sets the tuple <instance>.SOLVER
  <instance>.initproblem  -- sets the tuple <instance>.PROB
  <instance>.solve        -- solves the delay differential equation.
  <instance>.dde          -- does all of the above in one

In addition the module provides the functions:
    
  pydde.pastvalue    -- for accessing the history of the integration
  pydde.pastgradient -- for accessing the history of the derivative
  
After solving the DDE with <instance>.solve() the user may be interested in:
    
  <instance>.data -- contains the times (first column) and state variables
  
'''

import copy
import sys
import collections

try:
    from Numeric import array, zeros, concatenate
except ImportError:
    try:
        from numpy import array, zeros, concatenate
    except ImportError:
        print("Error: could not import either Numeric or SciPy.")  
        print("Please check that at least one of these is installed before using PyDDE.")
        raise ImportError('PyDDE.pydde')

try:
    import PyDDE.ddesolve as _ddesolve
except ImportError:
    print("Error: could not import the ddesolve integration routines!  The solver is non-functional.")
    _ddesolve = None
    raise ImportError('PyDDE.pydde')

    
class dde:

    def __init__(self):
        
        try:
            self._solved = 0
            self.data = None
            self.t = None
            self._set = [0,0]
            self.PROB = None
            self.SOLVER = None
        except:
            print("DDE Warning: something went wrong during initialisation.")
            raise

        
    
    def initproblem(self, 
                    no_vars=1, no_cons=1,
                    nlag=1, nsw=1,
                    t0=0, t1=1,
                    initstate=array([1.0]), c=array([1.0]), otimes=array([0.0,1.0]),
                    grad=(lambda s,c,t: array([-c[0]*s[0]])),
                    switchfunctions=(lambda s,c,t: array([1.0])),
                    maps=(lambda s,c,t,swno: (s,c)),
                    storehistory=(lambda g,s,c,t: (s,g))):
        '''
        Sets the parameters, initial conditions, derivative and switch functions,
        and a host of other values.
        
        The purpose of this function is to ensure that the problem is sufficiently
        well-specified for the solver to do its job.
        
        Some type checking and conversion is carried out, but don't rely on it.
        
        Arguments:
        
           self,
           no_vars         = 1,
           no_cons         = 1,
           nhv             = 0
           nlag            = 0,
           nsw             = 0,
           initstate       = array([1]),
           c               = array([1]),
           grad            = (lambda s,c,t: array([-c[0]*s[0]])),
           switchfunctions = (lambda s,c,t: array([1.0])),
           maps            = (lambda s,c,t,swno: array(s)),
           storehistory    = (lambda g,s,c,t: (array(s),array(g)))):
        
        '''
        
        # Check that the supplied functions are callable.
        try:
            if not(isinstance(grad, collections.Callable)):
                raise TypeError
            if not(isinstance(switchfunctions, collections.Callable)):
                raise TypeError
            if not(isinstance(maps, collections.Callable)):
                raise TypeError
            if not(isinstance(storehistory, collections.Callable)):
                raise TypeError
        except TypeError as errstr:
            print("DDE Error: User supplied function not callable:", errstr)
            print("DDE Error: Problem initialisation failed!")
            return 0
        
        # Check that the number of constants and variables are reasonable.
        try:
            if not(no_vars > 0):
                raise TypeError
        except TypeError as errstr:
            print("DDE Error: Number of state variables not positive:", errstr)
            return 0
        
        if (no_cons != len(c)):
            print("DDE Warning: Number of constants no_cons reset to agree with length of constants array c.")
            no_cons = len(c)
        
        # check that the state scale is of the appropriate length
        if (self.SOLVER != None):
            if (len(self.SOLVER[-1]) < no_vars):
                tempstsc = zeros((no_vars), 'd')
                for i in range(len(self.SOLVER[-1])):
                    tempstsc[i] = self.SOLVER[-1][i]
                self.SOLVER = (self.SOLVER[:-1],tempstsc)
            else:
                self.SOVLER = (self.SOLVER[:-1],self.SOLVER[-1][:no_vars])
        
        # Check that the outputs of the supplied functions have the right length,
        # if they will be used.
        try:
            g = grad(initstate,c,t0)
            if (len(g) != no_vars):
                raise IndexError
            if (nsw > 0):
                sw = switchfunctions(initstate,c,t0)
                if (len(sw) != nsw):
                    raise IndexError
                for i in range(nsw):
                    ms,mc = maps(initstate,c,t0,i)
                    if ((len(ms) != no_vars) or (len(mc) != no_cons)):
                        raise IndexError
        except IndexError as errstr:
            print("DDE Error: Output of user supplied function incorrect:", errstr)
            print("DDE Error: Problem initialisation failed!")
            return 0
        
        try:
            self.PROB = (int(no_vars), int(no_cons), # number of variables and constants
                         int(no_vars), int(nlag), int(nsw), int(len(otimes)),
                         float(t0), float(t1),
                         array(list(map(float,initstate))), array(list(map(float,c))), array(list(map(float,otimes))),
                         grad, switchfunctions, maps, storehistory)
        except:
            print("DDE Error: Something is wrong: perhaps one of the supplied variables has the wrong type?")
            print("DDE Error: Problem initialisation failed!")
            return 0
        
        # all went well
        self._set[0] = 1
        return 1
    
        
    def initsolver(self, tol=1e-8, hbsize=10000, 
                   dt=0.1,
                   statescale=array([0.0])):
        '''
        Sets the tolerance, buffer size, initial step size and a scaling
        parameter to help deal with states close to 0.
        
        Arguments:
           
           tol                       = 0.000005,
           hbsize                    = 1000,
           OBSOLETE: datasize        = 1000,
           OBSOLETE: dout            = 0.01,
           dt                        = 1.0,
           OBSOLETE: mindtx          = 10.0**(-9),
           OBSOLETE: maxdtx          = 100,
           statescale                = array([0]).
        '''
        
        # some simple checks
        if (hbsize == None): hbsize = 10000
        if (hbsize <= 0): hbsize = 1  #if 0 or less, C code crashes
        if (dt == None): dt = 0.1
        
        # check that the state scale is of the appropriate length
        if (self.PROB != None):
            if (len(statescale) < self.PROB[0]):
                tempstsc = zeros((self.PROB[0]), 'd')
                for i in range(len(statescale)):
                    tempstsc[i] = statescale[i]
                statescale = tempstsc
            else:
                statescale = statescale[:self.PROB[0]]
        
        try:
            self.SOLVER = (float(tol), int(hbsize),
                           float(dt), # output and integration timesteps
                           array(list(map(float,statescale))))
        except:
            print("DDE Error: Something is wrong: perhaps one of the supplied variables has the wrong type?")
            print("DDE Error: Solver initialisation failed!")
            return 0
        
        # all went well
        self._set[1] = 1
        return 1
        
    
    def solve(self, again=0):
        try:
            assert(self._set[0]*self._set[1])
            
            if (not(self._solved) or again):
                
                dd = _ddesolve.dde(self.PROB, self.SOLVER)
                self._solved = 1;
                
                # Attempt to get around a conflict arising when two dde() objects
                # are being manipulated.  This could be due to persistence of some
                # global variables in ddesolve but it's not clear.
                # Better solutions will be most welcome!
                self.data = copy.deepcopy(dd)
                del dd
                
                #clean(1)
                
            else:
                print(self._solved)
                print(again)
                print("DDE Error: Solver thinks the solution is already given in the <instance>.data.\n \
                       To force the solver to run again, use <instance>.solve(again=1)")
        
        except AssertionError:
            print("DDE Error: The DDE has not been properly initialised!")
            self._solved = 0
        except:
            print("DDE Error: Solution failed!")
            self._solved = 0
            raise
        
        return self.data
        
        
    def dde(self, 
            y=array([1]), 
            times=array([0.0,1.0]), 
            func=(lambda s,c,t: array([0.0])), 
            parms=array([1]), 
            switchfunc=(lambda s,c,t: array([1.0])),
            mapfunc=(lambda s,c,t,swno: (s,c)),
            tol=1e-8, 
            dt=0.1, 
            hbsize=10000,
            nlag=0, nsw=0, ssc=array([0.0])):
        '''
        Analogue of dde from ddesolve in R.  The last three arguments differ
        to give a little more control to the user.
        '''
        
        self.initproblem(no_vars=len(y), no_cons=len(parms),
                         nlag=nlag, nsw=nsw,
                         t0=times[0], t1=times[-1],
                         initstate=y, c=parms, otimes=times,
                         grad=func,
                         switchfunctions=switchfunc,
                         maps=mapfunc,
                         storehistory=(lambda g,s,c,t: (s,g)))
        self.initsolver(tol=tol, hbsize=hbsize, dt=dt, statescale=ssc)
        
        return self.solve(again=1)



def pastvalue(i,t,j):
    return _ddesolve.pastvalue(int(i),float(t),int(j))
    
def pastgradient(i,t,j):
    return _ddesolve.pastgradient(int(i),float(t),int(j))

def clean(wipe):
    try:
        assert(_ddesolve.clean(wipe))
    except AssertionError:
        print("DDE Error: Problem cleaning up after the solver.")
        return 0
    except:
        return 0
    return 1

