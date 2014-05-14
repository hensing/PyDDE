/***************************************************************************/
/* Original code by Simon N. Wood, 1999.                                   */
/* Ported from ddesolve for R by Alex Couture-Beil, Jon T. Schnute, and    */
/* Rowan Haigh, 2007.                                                      */
/* Updated for PyDDE and ddesolve by Benjamin J. Cairns, 2007-2008.        */
/* <ben.cairns@ceu.ox.ac.uk>                                               */
/***************************************************************************/

#ifndef _DDESOLVE95_PY_H_
#define _DDESOLVE95_PY_H_ 1

#include <stdio.h>

typedef struct { 
	int no_var, no_otherVars;
	int nhv,nlag,nsw;
	double dt,t0,t1,tol;
	long hbsize;
	char **cname,*initialtext,*initialtitle,**cinfo;
	FILE *file;
	int quit,newrun,cont,*findex,fileno;
	double **vals, *tmp_other_vals;
	int vals_size, vals_ind;
	double current_t;
	double *otimes; /* bjc 2007-05-08*/
	int no_otimes; /* bjc 2007-05-08*/
} globaldatatype;

/* globaldatatype data; MOVED to ddesolve95.c */
/* Thanks to Josh Lippai of Pomona College for this fix */

void numerics(double *c,int cont);
void setupglobaldata(int no_vars, int no_switch, double *settings, double *otimes, int no_otimes);
void freeglobaldata(void);

#endif /* _DDESOLVE95_PY_H_ */
