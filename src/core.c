/*
  # File Quor/src/core.c
  # 
  # Quor package for R (http://www.R-project.org)
  # Copyright (C) 2014 Adriano Polpo, Carlos A. de B. Pereira, Cassio P. de Campos.
  #
  #    This program is free software: you can redistribute it and/or modify
  #    it under the terms of the GNU General Public License as published by
  #    the Free Software Foundation, either version 3 of the License, or
  #    (at your option) any later version.
  #
  #    This program is distributed in the hope that it will be useful,
  #    but WITHOUT ANY WARRANTY; without even the implied warranty of
  #    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  #    GNU General Public License for more details.
  #
  #    You should have received a copy of the GNU General Public License
  #    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <R.h> 
#include <Rinternals.h> 
#include <Rmath.h>

#include <stdlib.h>

#define DEBUG(x)
#define Inf 1.e50

double getp(SEXP t, SEXP d, int x, int y) {
  return REAL(t)[(x) + INTEGER(d)[0] * (y)];
}
int geti(SEXP t, SEXP d, int x, int y) {
  return INTEGER(t)[(x) + INTEGER(d)[0] * (y)];
}

SEXP quorccore(SEXP data_, SEXP m_, SEXP lpr_, SEXP upr_, SEXP ngroups_, SEXP maxm_, SEXP nrows_, SEXP permmat, SEXP verbose_) {
  char vname[100];
  int i, jnext, j, k, ires;
  double v, vv;
  SEXP dimdata1;
  SEXP dimdata2;
  SEXP data1;
  SEXP data2;
  SEXP w, rl, nm;
  SEXP dimlpr = getAttrib(lpr_, R_DimSymbol);
  SEXP dimupr = getAttrib(upr_, R_DimSymbol);
  SEXP dimpermmat = getAttrib(permmat, R_DimSymbol);
  DEBUG(Rprintf("dimlpr=%d dimupr=%d dimpermmat=%d\n",INTEGER(dimlpr)[0],INTEGER(dimupr)[0],INTEGER(dimpermmat)[0]));
  double *m = REAL(m_);
  int perm;
  int nperms = INTEGER(dimpermmat)[0];
  int ngroups = INTEGER(ngroups_)[0];
  int maxm = INTEGER(maxm_)[0];
  int nrows = INTEGER(nrows_)[0];
  int verbose = INTEGER(verbose_)[0];
  int whoisi, whoisi1;
  int **Jl;
  double **D;
  double ***mpr;
  double d1, d2;

  if(INTEGER(dimlpr)[0]!=ngroups || INTEGER(dimupr)[0]!=ngroups) {
    REprintf("Wrong specification of lpr and upr - results will be wrong\n");
  }

  DEBUG(Rprintf("[debug] Number of orderings: %d\n",nperms);)
    DEBUG(Rprintf("[debug] Number of groups: %d\n",ngroups);)

    if(ngroups > 2) {
      //## When more than 2 groups, build also a cache of differences in cummulative binomials
      //## This takes in total O(ngroups*m^2) time, where m is the size of data for one covariate
      //## This step might be slow if we think of a single covariate, but it greatly speeds up
      //## computations in case there are many covariates (much more than ngroups) to be processed
      mpr = (double ***) malloc(sizeof(double **) * ngroups);
      for(i = 0; i < ngroups; i++) {
	mpr[i] = (double **) malloc(sizeof(double *) * (m[i]+1));
	for(j = 0; j < m[i]+1; j++) {
	  mpr[i][j] = (double *) malloc(sizeof(double) * m[i]);
	  for(k = 0; k < m[i]; k++) {
	    if(k <= j) mpr[i][j][k] = -Inf;
	    else {
	      d1 = getp(lpr_,dimlpr,i,k);
	      d2 = getp(lpr_,dimlpr,i,j);
	      if(d1 <= -0.1*Inf) mpr[i][j][k] = -Inf;
	      //## log( exp(lpr[i,k]) - exp(lpr[i,j]) )
	      else if (d2 <= -0.1*Inf) mpr[i][j][k] = d1;
	      else mpr[i][j][k] = d1 + log(1.0 - exp(d2 - d1));
	    }
	  }
	}
      }
    }

  Jl = (int **) malloc(sizeof(int *) * (ngroups-1));
  for(i = 0; i < ngroups-1; i++)
    Jl[i] = (int *) malloc(sizeof(int) * maxm);
  D = (double **) malloc(sizeof(double *) * ngroups);
  for(i = 0; i < ngroups; i++)
    D[i] = (double *) malloc(sizeof(double) * maxm);

  PROTECT(rl = allocVector(VECSXP,nperms));
  PROTECT(nm = allocVector(STRSXP,nperms));
  setAttrib(rl, R_NamesSymbol, nm);
  for(perm = 0; perm < nperms; perm++) {
    if(verbose) Rprintf("[info] Permutation %d started\n",perm+1);
    sprintf(vname, "result.%d", perm+1);
    PROTECT(w = allocVector(REALSXP, nrows));
    SET_VECTOR_ELT(rl, perm, w);
    UNPROTECT(1);
    SET_STRING_ELT(nm, perm, mkChar(vname));
    for(ires = 0; ires < nrows; ires++) {
      if(verbose && ires > 0 && ires % 100000 == 0) Rprintf("[info] %d covariates done so far\n",ires);

      //## First, for each element in a group i, mark who is the closest element
      //## to it in the next group i+1. The loop takes O(m) time in total. The
      //## allocation takes O(m*ngroups) space, which is linear in the input
      for(i = 0; i < ngroups-1; i++) {
	whoisi = geti(permmat,dimpermmat,perm,i)-1;
	whoisi1 = geti(permmat,dimpermmat,perm,i+1)-1;
	DEBUG(Rprintf(": whoisi=%d whoisi1=%d\n",whoisi,whoisi1));
	data1 = VECTOR_ELT(data_,whoisi);
	data2 = VECTOR_ELT(data_,whoisi1);
	dimdata1 = getAttrib(data1, R_DimSymbol);
	dimdata2 = getAttrib(data2, R_DimSymbol);
	if (INTEGER(dimdata1)[0]!=nrows || INTEGER(dimdata2)[0]!=nrows)
	  REprintf("Wrong data specification - results will be wrong\n");
	jnext = 0;
	for(j = 0; j < m[whoisi]; j++) {
	  while(jnext < m[whoisi1] && getp(data1,dimdata1,ires,j) >= getp(data2,dimdata2,ires,jnext)) jnext++;
	  Jl[i][j] = jnext;
	}
      }

      //## The dynamic programming begins. 
      //## The first line gets simply the prob of the quantile to the left of j, for each j
      whoisi = geti(permmat,dimpermmat,perm,0)-1;
      DEBUG(Rprintf("+ whoisi=%d\n",whoisi);)
	for(j = 0; j < m[whoisi]; j++) {
	  D[0][j] = getp(lpr_,dimlpr,whoisi,j);
	  DEBUG(Rprintf("D[%d][%d]=%.16lf\n",0,j,D[0][j]));
	}
      for (i = 1; i < ngroups; i++)
	for(j = 0; j < maxm; j++) D[i][j] = -Inf;

      //## If there only 2 groups, next loop is skipped. If run, it takes O(m^2) time in total
      for (i = 1; i < ngroups-1; i++) {
	whoisi = geti(permmat,dimpermmat,perm,i)-1;
	whoisi1 = geti(permmat,dimpermmat,perm,i-1)-1;
	DEBUG(Rprintf("@ whoisi=%d whoisi1=%d\n",whoisi,whoisi1);)
	  for (j = 0; j < m[whoisi]; j++) {
	    //## For each group, for each position j as separator, take the best solution for
	    //## all previous groups that are completely to the left of the separator, plus
	    //## the chance that the quantile of the current group will fall between that best
	    //## position for all previous groups and the separator j
	    v = D[i-1][0] + mpr[whoisi][Jl[i-1][0]][j];
	    DEBUG(k=0; Rprintf("$ vv[%d]=%.16lf D=%.16lf J=%d j=%d, mpr=%.16lf\n",k,vv,D[i-1][k],Jl[i-1][k],j,mpr[whoisi][Jl[i-1][k]][j]));
	    for(k = 1; k < m[whoisi1]; k++) {
	      vv = D[i-1][k] + mpr[whoisi][Jl[i-1][k]][j];
	      DEBUG(Rprintf("$ vv[%d]=%.16lf D=%.16lf J=%d j=%d, mpr=%.16lf\n",k,vv,D[i-1][k],Jl[i-1][k],j,mpr[whoisi][Jl[i-1][k]][j]));
	      if(vv > v) v = vv;
	    }
	    D[i][j] = v;
	    DEBUG(Rprintf("$ whoisi=%d whoisi1=%d perm=%d var=%d is %.16lf\n",whoisi,whoisi1,perm,ires,v));
	  }
      }

      //## Compose the result as the best solution found until the previous group with the prob
      //## that the final quantile falls to the right of the separator (the best separator is
      //## obtained in the maximization. The line takes O(m) time
      whoisi = geti(permmat,dimpermmat,perm,ngroups-1)-1;
      whoisi1 = geti(permmat,dimpermmat,perm,ngroups-2)-1;
      DEBUG(Rprintf("! whoisi=%d whoisi1=%d\n",whoisi,whoisi1));
      v = D[ngroups-2][0] + getp(upr_,dimupr,whoisi,Jl[ngroups-2][0]);
      DEBUG(Rprintf("vv[%d]=%.16lf D=%.16lf J=%d upr=%.16lf\n",0,v,D[ngroups-2][0],Jl[ngroups-2][0], getp(upr_,dimupr,whoisi,Jl[ngroups-2][0])));
      for(k = 1; k < m[whoisi1]; k++) {
	vv = D[ngroups-2][k] + getp(upr_,dimupr,whoisi,Jl[ngroups-2][k]);
	DEBUG(Rprintf("vv[%d]=%.16lf D=%.16lf J=%d upr=%.16lf\n",k,vv,D[ngroups-2][k],Jl[ngroups-2][k],getp(upr_,dimupr,whoisi,Jl[ngroups-2][k])));
	if(vv > v) v = vv;
      }
      DEBUG(Rprintf("result perm=%d var=%d is %.16lf\n",perm,ires,v));
      REAL(w)[ires] = v;
    }
  }
   
  if(ngroups > 2) {
    for(i = 0; i < ngroups; i++) {
      for(j = 0; j < m[i]+1; j++)
	free(mpr[i][j]);
      free(mpr[i]);
    }
    free(mpr);
  }
  for(i = 0; i < ngroups-1; i++) free(Jl[i]);
  free(Jl);
  for(i = 0; i < ngroups; i++) free(D[i]);
  free(D);
  UNPROTECT(2);
  return rl;
}

