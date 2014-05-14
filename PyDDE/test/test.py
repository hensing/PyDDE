'''
The following are commands for running this file.  Change the
argument of the second line to the directory where this test
file resides.  For example, if your file is in
/usr/lib/python2.4/site-packages/PyDDE/test, then use:

  os.chdir("/opt/lib/python2.4/site-packages/PyDDE/test")

Good luck!

import os
os.chdir("/Volumes/Working/projects/pydde/PyDDE/test")
execfile("test.py")

'''
    
try:
    from scipy import *
except ImportError:
    try:
        from numpy import *
    except ImportError:
        try:
            from Numeric import *
        except:
            print "Could not import Numeric, numpy and scipy.  Test cannot run."
            raise StandardError, "PyDDE test script."

try:
    import PyDDE.pydde as p
except ImportError:
    print "Could not import PyDDE.pydde.  Test cannot run."
    raise StandardError, "PyDDE test script."

try:
    import timing
    timeit = 1
except ImportError:
    print "Could not import timing module.  No timing of tests will occur."
    timeit = 0
    
'''
ODE example from Solv95 distribution.

This model is a Lotka-Volterra predator-prey system.
'''

if timeit:
    timing.start()

try:
    del ode_eg
except:
    print "DDE Test: Running for the first time..."

ode_eg = p.dde()

odegrad = (lambda s, c, t: array( [ c[2]*s[0]-c[3]*s[0]*s[1], c[0]*s[0]*s[1]-c[1]*s[1] ] ) )
    
odecons = array([0.005, 0.2, 1.0, 0.02, 100.0, 100.0])
odeist = array([odecons[4],odecons[5]])

ode_eg.initproblem(no_vars=2, no_cons=6, nlag=0, nsw=0,
                   t0=0.0, t1=300.0,
                   initstate=odeist, c=odecons, otimes=arange(0.0,300.0,1.0),
                   grad=odegrad)

odestsc = array([0,0])

ode_eg.initsolver(tol=1*10**(-8), hbsize=0,
                  dt=1.0, 
                  statescale=odestsc)

ode_eg.solve()
#print ode_eg.data



'''
DDE example from Solv95 distribution.

This model is a model for Nicholson's (1954) blowflies, as given by Gurney and Nisbet (1981)
'''

try:
    del dde_eg
except:
    print "DDE Test: Running for the first time..."

dde_eg = p.dde()

def ddegrad(s, c, t):
    alag = 0.0
    if (t>c[0]):
        alag = p.pastvalue(0,t-c[0],0)
    return array( [ c[2]*alag*exp(-alag/c[3])-c[1]*s[0] ] )
    
def ddesthist(g, s, c, t):
    return (s, g)

ddecons = array([12.0,0.25,10.0,300.0,100.0])
ddeist = array([ddecons[4]])
ddestsc = array([0])

dde_eg.dde(y=ddeist, times=arange(0.0, 300.0, 1.0), 
           func=ddegrad, parms=ddecons, 
           tol=0.000005, dt=1.0, hbsize=1000, nlag=1, ssc=ddestsc)

#print(dde_eg.data)

'''
SDDE (DDE with switches) example from Solv95 distribution.
'''

try:
    del sdde_eg
except:
    print "DDE Test: Running for the first time..."

sdde_eg = p.dde()

def sddegrad(s, c, t):
    g = array([0.0,0.0])
    g[0] = -c[0]*s[0]*s[1]
    g[1] = -c[5]*s[1]
    if (t>c[4]):
        g[1] = c[3]*p.pastvalue(0,t-c[4],0)*p.pastvalue(1,t-c[4],0)-c[5]*s[1]
    return g

def sddesw(s, c, t):
    sw = [0.0,0.0]
    sw[0]=sin(2*pi*t/c[1])         # add resource
    sw[1]=sin(2*pi*(t-c[4])/c[1])  # step nicely around a discontinuity
    return array(sw)
    
def sddemaps(s, c, t, swno):
    if (swno==0):
        return (array([s[0]+c[2],s[1]]),c)
    else:
        return (s,c)

sddecons = array([0.1, 10.0, 50.0, 0.05, 5.0, 0.02])
sddeist = array([0.0, 1.0])
sddestsc = 0*sddeist
                  
sdde_eg.dde(y=sddeist, times=arange(0.0, 300.0, 1.0), 
           func=sddegrad, parms=sddecons, 
           switchfunc=sddesw, mapfunc=sddemaps,
           tol=0.000005, dt=1.0, hbsize=1000, nlag=1, nsw=2, ssc=ddestsc)

#print sdde_eg.data





odetitles = ["prey", "predators"]
sddetitles = [r"\large resources", r"\large consumers"]
ddetitles = ["blowflies"]

if timeit:
    timing.finish()
    print str(timing.seconds())+"."+str(timing.milli())+" seconds"
    
    
try:
    from pyx import *

    '''
    Print the ODE model
    '''
    odepp = graph.graphxy(width=15, height=12,
                          x=graph.axis.linear(min=0,
                                              max=300,
                                              title=r"time $t$"),
                          y=graph.axis.linear(min=0,
                                              max=220,
                                              title=r"population sizes"),
                          key=graph.key.key(pos="tl"))
    attribs = [[style.linestyle.solid,style.linewidth.thick,color.rgb.blue],
               [style.linestyle.solid,style.linewidth.thick,color.rgb.red]]
    for i in range(len(odetitles)):
        thisCol = ode_eg.data[:,i+1]
        thisCol = reshape(thisCol,(thisCol.shape[0],1))
        thist = reshape(ode_eg.data[:,0],(ode_eg.data[:,0].shape[0],1))
        temp = concatenate((thist,thisCol),1)
        temp2 = list()
        for j in range(temp.shape[0]):
             temp2.append(tuple(temp[j,:]))
        if odetitles[i] != '':
            odepp.plot(graph.data.list(temp2, x=1, y=2, title=r""+str(odetitles[i])),
                       [graph.style.line(lineattrs=attribs[i])])
    
    odepp.writeEPSfile("../test/ode_eg.eps")

    '''
    Print the DDE model
    '''
    ddepp = graph.graphxy(width=15, height=12,
                          x=graph.axis.linear(min=0,
                                              max=300,
                                              title=r"time $t$"),
                          y=graph.axis.linear(min=0,
                                              max=4000,
                                              title=r"population size"),
                          key=graph.key.key(pos="tl"))
    attribs = [[style.linestyle.solid,style.linewidth.thick,color.rgb.blue],
               [style.linestyle.dashed,style.linewidth.thick,color.rgb.red]]
    for i in range(len(ddetitles)):
        thisCol = dde_eg.data[:,i+1]
        thisCol = reshape(thisCol,(thisCol.shape[0],1))
        thist = reshape(dde_eg.data[:,0],(dde_eg.data[:,0].shape[0],1))
        temp = concatenate((thist,thisCol),1)
        temp2 = list()
        for j in range(temp.shape[0]):
            temp2.append(tuple(temp[j,:]))
        if ddetitles[i] != '':
            ddepp.plot(graph.data.list(temp2, x=1, y=2, title=r""+str(ddetitles[i])),
                       [graph.style.line(lineattrs=attribs[i])])

    ddepp.writeEPSfile("../test/dde_eg.eps")

    '''
    Print the SDDE model
    '''
    sddepp = graph.graphxy(width=15, height=12,
                           x=graph.axis.linear(min=0,
                                               max=300,
                                               title=r"\large time $t$"),
                           y=graph.axis.linear(min=0,
                                               max=200,
                                               title=r"\large states"),
                           key=graph.key.key(pos="tl"))
    attribs = [[style.linestyle.solid,style.linewidth.thick,color.rgb.blue],
               [style.linestyle.solid,style.linewidth.thick,color.rgb.red]]
    for i in range(len(odetitles)):
        thisCol = sdde_eg.data[:,i+1]
        thisCol = reshape(thisCol,(thisCol.shape[0],1))
        thist = reshape(sdde_eg.data[:,0],(sdde_eg.data[:,0].shape[0],1))
        temp = concatenate((thist,thisCol),1)
        temp2 = list()
        for j in range(temp.shape[0]):
            temp2.append(tuple(temp[j,:]))
        if sddetitles[i] != '':
            sddepp.plot(graph.data.list(temp2, x=1, y=2, title=r""+str(sddetitles[i])),
                        [graph.style.line(lineattrs=attribs[i])])
                
    sddepp.writeEPSfile("../test/sdde_eg.eps")
        
except ImportError:
    print "DDE Test: Could not import PyX!  Cannot print output from examples."
except:
    print "DDE Test: An error occured during attempt to plot data."
    raise
