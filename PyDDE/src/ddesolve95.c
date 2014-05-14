/***************************************************************************/
/* Original code by Simon N. Wood, 1999.                                   */
/* Ported from ddesolve for R by Alex Couture-Beil, Jon T. Schnute, and    */
/* Rowan Haigh, 2007.                                                      */
/* Updated for PyDDE and ddesolve by Benjamin J. Cairns, 2007-2008.        */
/* <ben.cairns@ceu.ox.ac.uk>                                               */
/***************************************************************************/

#include <stdlib.h>
#include "ddesolve95.h"
#include "ddeq.h"

#define CH_BUF_SIZE 128

globaldatatype data; /* MOVED from ddesolve95.h */
/* Thanks to Josh Lippai of Pomona College for this fix */

int the_test_phase=0;

/*===========================================================================*/
void info(char *str)
{
	wARNING(str);
}

/*===========================================================================*/
void output(double *s,double t)
{
	/*data.vals[0] is for t,
	  and [1..(no_var+1)] are reserved for s[0..no_var] vars
	  */
	int i;
	static double *dummy_var=NULL;
	if( dummy_var == NULL )
		dummy_var = malloc( data.no_var*sizeof(double) );
	data.vals[0][data.vals_ind] = t;
	for( i = 0; i < data.no_var; i++ )
		data.vals[i+1][data.vals_ind] = s[i];
	
	/*ACB hack - call grad to pull out any other data*/
	if( data.no_otherVars > 0 )
		grad(dummy_var,s,NULL,t);
	
	for( i = 0; i < data.no_otherVars; i++ )
		data.vals[1+data.no_var+i][data.vals_ind] = data.tmp_other_vals[i];
	
	data.vals_ind++;
	
	if (data.vals_ind >= data.vals_size) {
		for(i=0;i<(1+data.no_var+data.no_otherVars);i++) {
			data.vals[i] = (double*)realloc(data.vals[i], sizeof(double)*2*data.vals_size);
			if (data.vals[i]==NULL)
				eRROR("memory (re)allocation failed");
		}
		data.vals_size *= 2;
	}
}

/*===========================================================================*/
/* cont is zero for fresh run, and 1 for continuation */
void numerics(double *c,int cont)
{ 
	static double *s;
	double t0, t1, dt, *otimes; /* bjc 2007-05-08*/
	int ns, nsw, nhv, nlag, reset=1, fixstep=0, no_otimes; /* bjc 2007-05-08*/
	static int first=1;
	long hbsize;
	ns=data.no_var;
	nsw=data.nsw;
	nhv=data.nhv;
	nlag=data.nlag;
	t0=data.t0;
	t1=data.t1;
	dt=data.dt;
	hbsize=data.hbsize;
	otimes=data.otimes;
	no_otimes=data.no_otimes; /* bjc 2007-05-08*/
  
	if (cont) {
		reset=0;
	} else {
		if (!first) {
			free(s);
			first=0;
  		}
		s=(double *)calloc(data.no_var,sizeof(double));
		ddeinitstate(s,c,t0);
	}
	dde(s,c,t0,t1,&dt,data.tol,otimes,no_otimes,ns,nsw,nhv,hbsize,nlag,reset,fixstep); /* bjc 2007-05-08*/
	data.dt=dt;
}

/*===========================================================================*/
void setupglobaldata(int no_vars, int no_switch, double *settings, double *otimes, int no_otimes) /* bjc 2007-05-08*/
{ 
	int i;
	const int no_otherVars = 0;

	data.tol=settings[0];
	data.t0=settings[1];        /* start time */
	data.t1=settings[2];        /* stop time */
  
	data.dt=settings[3];        /* initial timestep */
	
	data.hbsize=settings[4];    /* how many past values to store for each history variable */
	data.no_var=no_vars;
  
	// OBSOLETE: data.no_otherVars=no_otherVars;

	data.nsw=no_switch;          /* number of switch varaibles */  
	data.nhv=no_vars;         /* Number of history (lagged) variables */
	data.nlag=1;        /* Number of lag markers per history variable (set to 1 if unsure)*/

	/* enter out times into the data structure */
	data.otimes = otimes; /* bjc 2007-05-08: could be NULL*/
	data.no_otimes = no_otimes; /* bjc 2007-05-08: >= 0  */

	data.vals_size=1000; /* size will grow, this is just initial min size */
	data.vals_ind=0;
	data.vals = (double**)malloc(sizeof(double*)*(data.no_var+no_otherVars+1));
	if (data.vals==NULL)
		eRROR("memory allocation failed");
	for(i=0;i<(data.no_var+no_otherVars+1);i++) {
		data.vals[i]=(double*)malloc(sizeof(double)*data.vals_size);
		if (data.vals[i]==NULL)
			eRROR("memory allocation failed");
	}
	if (data.no_otherVars>0) {
		data.tmp_other_vals = (double*)malloc(sizeof(double)*data.no_otherVars);
		if (data.tmp_other_vals==NULL) {
			eRROR("memory allocation failed");
		}
	} else {
		data.tmp_other_vals=NULL;
	}
}

/*===========================================================================*/
void freeglobaldata()
{
	int i;
	if (data.vals) {
		for(i=0;i<(data.no_var+data.no_otherVars+1);i++) {
			free(data.vals[i]);
		}
		free(data.vals);
		data.vals=NULL;
	}
	if (data.tmp_other_vals) {
		free(data.tmp_other_vals);
		data.tmp_other_vals=NULL;
	}
}

