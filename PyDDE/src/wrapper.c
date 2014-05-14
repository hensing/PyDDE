/***************************************************************************/
/* Some original code by Simon N. Wood, 1999.                              */
/* Updated for Python port by Benjamin J. Cairns, 2005-2008.               */
/* <ben.cairns@ceu.ox.ac.uk>                                               */
/***************************************************************************/

#include <Python.h>
#include <math.h>
#include <assert.h>
#include "ddeq.h"
#include "ddesolve95.h"

#ifndef PY_ARRAY_UNIQUE_SYMBOL
#define PY_ARRAY_UNIQUE_SYMBOL Py_Array_DDE
#endif

#include <numpy/arrayobject.h>

extern globaldatatype data; /* ddesolve95.c */

/***********************************************************/
/* Stuff for initial states and statescale (error control) */
/***********************************************************/

// Some variables it is convenient to make external to functions
static double *initstateDbl, *statescaleDbl, *cDbl, *otimesDbl;
static int no_cons, no_vars, nhv, nlag, nsw, no_otimes;
static long hbsize;
static double dt, t0, t1, tol;
extern int failed;

// The Python functions passed in by the user
static PyObject *grad_func = NULL;
static PyObject *switchfunctions_func = NULL;
static PyObject *map_func = NULL;
static PyObject *storehistory_func = NULL;

// Functions providing user-supplied data.
void ddeinitstate(double *s, double *c, double t)
/* initialise state variables and any global constants here, you can use c */
{
	int i;
	for(i=0;i<data.no_var;i++)
		s[i]=initstateDbl[i];
}

void statescale(double* scale) {
    scale = statescaleDbl;
}



/***************************************************************************/
/* Utility to convert ONE-DIMENSIONAL Python float arrays to double arrays */
/***************************************************************************/

double *dblArray_from_PyArray(PyObject *arrayPy) {
    PyArrayObject *contig_array;
    int n, i;
    double *array;
        
    // Next line handles scalars.and non-contiguous arrays
    contig_array = (PyArrayObject *)PyArray_ContiguousFromObject(arrayPy, PyArray_DOUBLE, 1, 0);
    //contig_array = (PyArrayObject *)arrayPy;
    if (contig_array == NULL) {
        return NULL;
    }
    else if (contig_array->nd > 1) {
        PyErr_SetString(PyExc_TypeError, "Array is not (at most) 1-dimensional!");
        Py_DECREF(contig_array);
        return NULL;
    }
    else {
        
        n = contig_array->dimensions[0];
        //printf("Allocating in dblArray_from_PyArray with %d elements.\n",  n);
        array = (double *)calloc(n, sizeof(double));
        for (i=0; i < n; i++) {
            array[i] = *(double *)(contig_array->data + i*contig_array->strides[0]);
        }
    }
    
    // Free up the reference to arrayPy or its copy
    Py_DECREF(contig_array);
    
    // Note receiver must free
    return array;
}

PyObject *pyArray_from_DblArray(double *dblArray, int size) {
    
    PyObject *pyArray;
    PyArrayObject *pyA;
    int i;
    
    pyArray = PyArray_FromDims(1,&size,PyArray_DOUBLE);
    //pyA = (PyArrayObject *)PyArray_ContiguousFromObject(pyArray, PyArray_DOUBLE, 1, 1);
    pyA = (PyArrayObject *)pyArray;
    
    for (i=0;i<size;i++) {
        *(double *)(pyA->data + \
                    i*pyA->strides[0]) = dblArray[i];
    }
    
    // Note receiver must also Py_DECREF
    return pyArray;
}

PyObject *pyArray2_from_DblArray2(double **dblArray, int *size) {
    
    PyObject *pyArray;
    PyArrayObject *pyA;
    int i,j;
    
    pyArray = PyArray_FromDims(2,size,PyArray_DOUBLE);
    //pyA = (PyArrayObject *)PyArray_ContiguousFromObject(pyArray, PyArray_DOUBLE, 2, 2);
    pyA = (PyArrayObject *)pyArray;
    
    for (j=0;j<size[0];j++) {
        for (i=0;i<size[1];i++) {
            *(double *)(pyA->data + \
                        i*pyA->strides[1] + \
                        j*pyA->strides[0]) = dblArray[i][j];
        }
    }
    
    // Note receiver must also Py_DECREF
    return pyArray;
}

// Utility function that was useful in debugging the above.
void print_PyArray(PyObject *array) {
    PyArrayObject *a;
    int i;
    
    a = (PyArrayObject *)array;
    
    //printf("-- This is print_PyArray.\n");
    //printf("-- ... nd = %d\n", a->nd);
    //printf("-- ... dim. 1 size = %d\n", a->dimensions[0]);
    //printf("-- ... [ ");
    for (i=0;i<a->dimensions[0];i++) {
        //printf("%f ",*(double *)(a->data + i*a->strides[0]));
    }
    //printf("]\n");
}



/****************************************************************/
/* These function provides access to pastgradient and pastvalue */
/****************************************************************/

// prototypes for external functions
extern double pastgradient(int i, double t, int markno);
extern double pastvalue(int i, double t, int markno);

// Wrapper for pastgradient
PyObject *wrap_pastgradient(PyObject *self, PyObject *args) {
    int j,k;
    double tt, result;
    if (!PyArg_ParseTuple(args,"idi",&j,&tt,&k)) {
        PyErr_SetString(PyExc_TypeError, "Could not parse arguments in 'pastgradient'!");
        return Py_BuildValue("f",0.0);
    }
    result = pastgradient(j, tt, k);
    return Py_BuildValue("f",result);
}

// Wrapper for pastvalue
PyObject *wrap_pastvalue(PyObject *self, PyObject *args) {
    int j,k;
    double tt, result;
    if (!PyArg_ParseTuple(args,"idi",&j,&tt,&k)) {
        PyErr_SetString(PyExc_TypeError, "Could not parse arguments in 'pastvalue'!");
        return Py_BuildValue("f",0.0);
    }
    result = pastvalue(j, tt, k);
    return Py_BuildValue("f",result);
}



/************************************/
/* Cleanliness is next to dogliness */
/************************************/

/* vextern void freeglobaldata(); */

// Wrapper for freeglobaldata
PyObject *wrap_freeglobaldata(PyObject *self, PyObject *args) {
    int wipeit;
    
    if (!PyArg_ParseTuple(args,"i",&wipeit)) {
        PyErr_SetString(PyExc_TypeError, "Could not parse arguments in 'clean'!");
        return Py_BuildValue("i",0);
    }
    
    //printf("Cleaning global data.\n");
    if (wipeit) {
        freeglobaldata();
    }
    return Py_BuildValue("i",1);
}



/*********************/
/* The main business */
/*********************/

// The main function for solving things
PyObject *wrap_dde(PyObject *self, PyObject *args) {
    // Variables from PROB
    PyObject *initstateObj, *statescaleObj, *cObj, *otimesObj;
    PyObject *gradObj, *switchfunctionsObj, *mapObj, *storehistoryObj;
    // Return variables
    int outdims[2];
    PyObject *dataPy=NULL;
    double settings[5];

    // Ensure variables initialised!
    failed = 0;
    
    //printf("This is dde\n");
        
    // get the arguments, two tuples giving problem and solver parameters
    if (!PyArg_ParseTuple(args,"(iiiiiiddOOOOOOO)(dldO)",\
                          &no_vars,&no_cons,&nhv,&nlag,&nsw,&no_otimes,\
                          &t0,&t1,\
                          &initstateObj,&cObj,&otimesObj,\
                          &gradObj,&switchfunctionsObj,&mapObj,&storehistoryObj,\
                          &tol,&hbsize,&dt,&statescaleObj)) {
        PyErr_SetString(PyExc_TypeError, "Could not parse arguments in 'dde'!");
        failed = 1;
    }
    statescaleDbl = dblArray_from_PyArray(statescaleObj);
    cDbl = dblArray_from_PyArray(cObj);
    initstateDbl = dblArray_from_PyArray(initstateObj);
    otimesDbl = dblArray_from_PyArray(otimesObj);
    
    //printf("no_vars = %d, no_cons = %d, nhv = %d\n", no_vars, no_cons, nhv);
    //printf("dt = %f\n", dt);
    
    // Check that the user-supplied functions are callable and set them up for callbacks
    if (!PyCallable_Check(gradObj)) { // gradient
        PyErr_SetString(PyExc_TypeError, "The 'grad' function must be a callable object!");
        failed = 1;
    }
    else {
        grad_func = gradObj; Py_INCREF(grad_func); 
        //printf("grad_func = %p\n",grad_func);
    }
    if (!PyCallable_Check(switchfunctionsObj)) { // switchfunctions
        PyErr_SetString(PyExc_TypeError, "The 'switchfunctions' function must be a callable object!");
        failed = 1;
    }
    else {
        switchfunctions_func = switchfunctionsObj; Py_INCREF(switchfunctions_func); 
        //printf("switchfunctions_func = %p\n",switchfunctions_func);
    }
    if (!PyCallable_Check(mapObj)) { // map
        PyErr_SetString(PyExc_TypeError, "The 'maps' function must be a callable object!");
        failed = 1;
    }
    else {
        map_func = mapObj; Py_INCREF(map_func); 
        //printf("map_func = %p\n",map_func);
    }
    if (!PyCallable_Check(storehistoryObj)) { // storehistory
        PyErr_SetString(PyExc_TypeError, "The 'storehistory' function must be a callable object!");
        failed = 1;
    }
    else {
        storehistory_func = storehistoryObj; Py_INCREF(storehistory_func); 
        //printf("storehistory_func = %p\n",storehistory_func);
    }
    
    // Initialise settings
    settings[0] = tol;     //printf("tol = %f\n", tol);
    settings[1] = t0;      //printf("t0 = %f\n", t0);
    settings[2] = t1;      //printf("t1 = %f\n", t1);
    settings[3] = dt;      //printf("dt = %f\n", dt);
    settings[4] = hbsize;  //printf("hbsize = %d\n", hbsize);
    
    // Now to the meat of the function...
    //printf("Data successfully loaded.\n");
    setupglobaldata(no_vars, nsw, settings, otimesDbl, no_otimes);
    numerics(cDbl, 0);
    if (failed) mESSAGE("Integration failed!  Outputs are None.");
    else mESSAGE("Integration completed without error.");
    
    if (failed) { // Something went wrong (user should see why) so return None.
        Py_INCREF(Py_None); // Be       nice to the reference count.
    //    Py_INCREF(Py_None); // Be extra nice to the reference count.
        dataPy = Py_None;
    //    tPy = Py_None;
    }
    else { // Nothing went wrong.  Do whatever necessary to wrap it up.
        
        outdims[0] = data.vals_ind;
        outdims[1] = no_vars+1;
        //printf("... data.outcount = %d\n", data->outcount);
        
        dataPy = pyArray2_from_DblArray2(data.vals, outdims); //printf("dataPy = %p, ", dataPy);
        //OBSOLETE: tPy = pyArray_from_DblArray(data->t, data->outcount); //printf("tPy = %p\n", tPy);
        
    }
    
    freeglobaldata();
    
    // Free memory.  DO NOT DECREMENT OBJECT REFERENCES FROM PyArgs_ParseTuple()
    free(initstateDbl); free(statescaleDbl); free(cDbl); free(otimesDbl);
    /* DO NOT wrap_freeglobaldata(1); this is done on the Python side by calling clean(1) */
    Py_DECREF(grad_func);Py_DECREF(switchfunctions_func);Py_DECREF(map_func);Py_DECREF(storehistory_func);
    
    return Py_BuildValue("N",dataPy);
}



/************************************************************/
/* Critical user-supplied functions and their mappings to C */
/************************************************************/

// Read the comments in the code of this function to understand the following 3
void switchfunctions(sw,s,c,t)
double *sw,*s,*c,t;
/* This routine sets the values of the switch functions. When the switch
	functions pass through zero from positive to negative the state variables
	may be reset in function map(). The switch functions should pass smoothly
	through 0 and should never have both value and first derivative zero. The
	same switch must not pass through zero from positive to negative more than
	once in an integration timestep. An example of a switch function is:
						sw[0]=sin(pi*t/30.0)
	which passes through zero every 60 time units. Switches may include state
	variables provided the above conditions are met. Note that to get 'Solver'
	style switches define twice as many switches and let e.g. sw[1]=-sw[0] */
{ 
    PyObject *result,*arglist;
    PyObject *sPy,*cPy;
    double *currentdata;
    int i;
    
    //printf("This is switchfunctions.\n");
    
    // Get the double arrays into Python arrays
    sPy = pyArray_from_DblArray(s,no_vars);
    cPy = pyArray_from_DblArray(c,no_cons);
    
    // Construct the argument list to pass back to the user-supplied function
    arglist = Py_BuildValue("(OOf)",sPy,cPy,t);
    
    // Obtain the result, convert back to a double array and put it in the right place.
    result = PyEval_CallObject(switchfunctions_func,arglist);
    currentdata = dblArray_from_PyArray(result); // put the results somewhere
    for (i=0;i<nsw;i++) {
        sw[i] = currentdata[i]; // copy across to sw
    }

    // Look after the reference counts to local variables
    Py_XDECREF(result);
    Py_DECREF(arglist);
    Py_DECREF(sPy);
    Py_DECREF(cPy);
    //printf("Freeing in switchfunctions.\n");
    free(currentdata);
}

void map(s,c,t,swno)
double *s,*c,t;int swno;
/* This routine is called whenever one of the switch functions passes through
	zero. 'swno' is the number of the switch function. The state variables
	can be changed discontinuously within this routine. eg:
   if (swno==1)
	  { s[0]=coeff[1]*(s[0]);}
	time and the coefficients should not be changed.
*/
{ 
    PyObject *result,*arglist;
    PyObject *sPy,*cPy,*snewPy,*cnewPy;
    double *currentdata;
    int i;
    
    //printf("This is map.\n");
    
    sPy = pyArray_from_DblArray(s,no_vars);
    cPy = pyArray_from_DblArray(c,no_cons);
    
    arglist = Py_BuildValue("(OOfi)",sPy,cPy,t,swno);
    
    result = PyEval_CallObject(map_func,arglist);
    //printf("result = %p\n", result);
    if (!PyArg_ParseTuple(result,"OO",&snewPy,&cnewPy)) {
        PyErr_SetString(PyExc_TypeError, "Could not parse results of 'map'!");
        Py_XDECREF(result);
        Py_DECREF(sPy);
        Py_DECREF(cPy);
        Py_DECREF(arglist);
        return;
    }
    else {
        
        currentdata = dblArray_from_PyArray(snewPy);
        for (i=0; i < no_vars; i++) {
            s[i] = currentdata[i];
        }
        
        free(currentdata);
        currentdata = dblArray_from_PyArray(cnewPy);
        for (i=0; i < no_cons; i++) {
            c[i] = currentdata[i];
        }

        // Look after the reference counts to local variables
        Py_XDECREF(result);
        Py_DECREF(arglist);
        Py_DECREF(sPy);
        Py_DECREF(cPy);
        //printf("Freeing in map.\n");
        free(currentdata);
    }
}

void grad(g,s,c,t)
double *g,*s,*c,t;
/* This routine must provide the gradients g for the state variables s.
	So ds[i]/dt=g[i]=fi(s,c,t) where c is the coefficient vector. lagged
   variables may be accessed here using pastvalue(i,x,j) which returns the
   ith (starting at zero) lagged variable at time x, using lag pointer k

   (lag pointers are used by pastvalue to store the history buffer location
    corresponding to a lag in order to save exectution time. For example if
    your code requires lagged varaible 0 at lags T and 2T for each gradient
    calculation then it is efficient to obtain these values using:
    pastvalue(0,t-T,0) and pastvalue(0,t-2*T,1) rather than
    pastvalue(0,t-T,0) and pastvalue(0,t-2*T,0). The latter works, it's just
    slower because more time is spent searching for lagged values)
*/
{ 
    PyObject *result,*arglist;
    PyObject *sPy,*cPy;
    double *currentdata;
    int i;
    
    //printf("This is grad.\n");
    
    sPy = pyArray_from_DblArray(s,no_vars);
    cPy = pyArray_from_DblArray(c,no_cons);
    
    arglist = Py_BuildValue("(OOf)",sPy,cPy,t);
    //printf("Got argument list.\n");    
    
    //printf("grad_func = %p\n", grad_func);
    //printf("arglist = %p\n", arglist);
    //PyObject_Print(grad_func, stdout, 0);
    //PyObject_Print(arglist, stdout, 0);
    //printf("Continuing...\n");
    
    //printf("grad_func = %p\n", grad_func);
    //printf("arglist = %p\n", arglist);
    result = PyEval_CallObject(grad_func,arglist);
    //printf("Got grad_func result.\n");
    assert(result);
    assert(PyArray_Check(result));
    //printf("result = %p\n", result);
    currentdata = dblArray_from_PyArray(result);
    //printf("Ready to copy gradient.\n");
    for (i=0;i < no_vars;i++) {
        g[i] = currentdata[i];
    }
        
    // Look after the reference counts to local variables
    Py_XDECREF(result);
    Py_DECREF(arglist);
    Py_DECREF(sPy);
    Py_DECREF(cPy);
    //printf("Freeing in grad.\n");
    free(currentdata);
}

void storehistory(his,ghis,g,s,c,t)
double *his,*ghis,*g,*s,*c,t;
/* This is the routine in which the values of the history variables at time
	t are calculated and put in the array his, along with gradients in ghis,
	using state variables s, gradients of s, g, and coefficients c
   e.g. if the state variable 2 is history variable 0, you would need the line:
   his[0]=s[2];ghis[0]=g[2];
*/
{ 
    PyObject *result,*arglist;
    PyObject *sPy,*cPy,*gPy;
    PyObject *hisPy = NULL,*ghisPy = NULL;
    double *currentdata;
    int i;
    
    //printf("This is storehistory.\n");
    
    sPy = pyArray_from_DblArray(s,no_vars);
    cPy = pyArray_from_DblArray(c,no_cons);
    gPy = pyArray_from_DblArray(g,no_vars);
    
    arglist = Py_BuildValue("(OOOf)",gPy,sPy,cPy,t);
    
    result = PyEval_CallObject(storehistory_func,arglist);
    if (!PyArg_ParseTuple(result,"OO",&hisPy,&ghisPy)) {
        PyErr_SetString(PyExc_TypeError, "Could not parse results of 'storehistory'!");
        Py_XDECREF(result);
        Py_DECREF(sPy);
        Py_DECREF(cPy);
        Py_DECREF(gPy);
        Py_DECREF(arglist);
        return;
    }
    else {
        
        currentdata = dblArray_from_PyArray(hisPy);
        for (i=0;i < nhv;i++) {
            his[i] = currentdata[i];
        }
        //printf("Freeing in storehistory.\n");
        free(currentdata);
        
        currentdata = dblArray_from_PyArray(ghisPy);
        for (i=0;i < nhv;i++) {
            ghis[i] = currentdata[i];
        }
        
        //printf("Look after the reference counts to local variables\n");
        Py_XDECREF(result);
        Py_DECREF(arglist);
        Py_DECREF(sPy);
        Py_DECREF(cPy);
        Py_DECREF(gPy);
        //printf("Freeing in storehistory.\n");
        free(currentdata);
        //printf("this was storehistory.\n");
    }
}



/****************************************/
/* All the required Python module stuff */
/****************************************/

// method definitions
static PyMethodDef ddesolveMethods[] = { 
         { "pastgradient", wrap_pastgradient, METH_VARARGS },
         { "pastvalue", wrap_pastvalue, METH_VARARGS },
         { "clean", wrap_freeglobaldata, METH_VARARGS },
         { "dde", wrap_dde, METH_VARARGS },
         { NULL, NULL } 
}; 

// module initialisation
void initddesolve(void) { 
         PyObject *m; 
         m = Py_InitModule("ddesolve", ddesolveMethods); 
         import_array();
} 

