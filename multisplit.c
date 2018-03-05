/* multisplit */
/* Author: F Tilmann */
/* Contact: tilmann|a|gfz-potsdam.de */

/* Create splitting estimate and error surface by grid search for splitting parameters */
/* Several different methods of measuring splitting are implemented */

/* (C) 2004 F Tilmann */
/* This source code is released under the GNU public license */

/* Code uses sac,  gsl and gslblas libraries and needs GMT programs to be installed */

/* History: */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stddef.h>
#include <math.h>

#include <gsl/gsl_vector_float.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_fft_real_float.h>
#include <gsl/gsl_fft_halfcomplex_float.h>
#include <gsl/gsl_blas.h>
/* #include <gsl/gsl_cblas.h> */
#include <gsl/gsl_eigen.h>

/* #include "saclib.h" */
/* multisplit.h includes: sac.h */
#include "sac_help.h"
#include "multisplit.h"  



int verbose=0;
#define VRB(command) { if(verbose) { command  ; fflush(stdout); }}

char warn_str[1024];
char abort_str[1024];

int main(int argc, char **argv)
{
  ms_params *par=(ms_params *) malloc(sizeof(ms_params));
  /* window definition */
  long wbeg,wlen;
  long refwbeg,refwlen;
  float maxsumlag;
  gsl_vector_float *hor1,*hor2,*ref1,*ref2;
 
  hor_split *hsplit;
  /* Grid search arrays */
  gsl_matrix *mat_res,*mat_pol,*mat_delay, *mat_alpha;

  long ldum1,ldum2;
  int n,m;
  int mod_par;
  char tmpstring[256],phase[9];

/*    invfisher test  */
/*   printf("1. %f\n",invfisher(2,20.85,0.68)); */
/*   printf("2. %f\n",invfisher(20.85,2,0.68)); */
/*   printf("3. %f\n",invfisher(2,20.85,0.95)); */
/*   printf("4. %f\n",invfisher(20.85,2,0.95)); */
/*   printf("PER(360,0)=%d\n",PER(360,0)); */
/*   printf("PER(360,340)=%d\n",PER(360,340)); */
/*   printf("PER(360,450)=%d\n",PER(360,450)); */
/*   printf("PER(360,-20)=%d\n",PER(360,-20)); */
/*   printf("PER(360,-380)=%d\n",PER(360,-380)); */
/*   exit(0); */

  parse(argc,argv, par);
  /* Check that sac files are consistent */
  if(check_consistency(par->hdr_hor1,par->hdr_hor2, 
                       CONSISTENCY_STATION | CONSISTENCY_EVENT | CONSISTENCY_VERBOSE | CONSISTENCY_BEGIN )) 
    abort_msg("Data files inconsistent");     
  make_rhs(&par->hdr_hor1,&par->data_hor1,&par->hdr_hor2,&par->data_hor2);

  if(!F_EQ(par->hdr_hor1->cmpaz,0.0)){ 
    printf("Rotating components from %f/%f to NE system\n",par->hdr_hor1->cmpaz,par->hdr_hor2->cmpaz);
    gsl_float_rotate(par->data_hor1,par->data_hor2,par->hdr_hor1->cmpaz);
    par->hdr_hor1->cmpaz=0.0;
    par->hdr_hor2->cmpaz=90.0;
  }

  if(par->method==CORREL) {
    if(check_consistency(par->hdr_hor1,par->method_q.cor_par.hdr_ref1,
			 CONSISTENCY_EVENT | CONSISTENCY_VERBOSE ))
      abort_msg("Inconsistency between data and reference files");
    if(check_consistency(par->method_q.cor_par.hdr_ref1,par->method_q.cor_par.hdr_ref2, 
		      CONSISTENCY_STATION | CONSISTENCY_EVENT | CONSISTENCY_VERBOSE | CONSISTENCY_BEGIN ))
      abort_msg("Reference files inconsistent");
    make_rhs(&par->method_q.cor_par.hdr_ref1,&par->method_q.cor_par.data_ref1,
	     &par->method_q.cor_par.hdr_ref2,&par->method_q.cor_par.data_ref2);
    if(!F_EQ(par->method_q.cor_par.hdr_ref1->cmpaz,0.0)){ 
      printf("Rotating components from %f/%f to NE system\n",par->method_q.cor_par.hdr_ref1->cmpaz,par->method_q.cor_par.hdr_ref2->cmpaz);
      gsl_float_rotate(par->method_q.cor_par.data_ref1,par->method_q.cor_par.data_ref2,par->method_q.cor_par.hdr_ref1->cmpaz);
      par->method_q.cor_par.hdr_ref1->cmpaz=0.0;
      par->method_q.cor_par.hdr_ref2->cmpaz=90.0;
    }
  }

  /* maxsumlag: maximum shift between seismograms in s */
  maxsumlag=par->model_q.split_par.top.timemax + par->model_q.split_par.bot.timemax ;

  if(par->method==CORREL) {
    maxsumlag += par->method_q.cor_par.maxshift;
  }
  hor1=find_window(par->hdr_hor1,par->data_hor1,&par->window,maxsumlag,&wbeg,&wlen,phase);
  hor2=find_window(par->hdr_hor2,par->data_hor2,&par->window,maxsumlag,&ldum1,&ldum2,tmpstring);
  if (ldum1 != wbeg || ldum2 != wlen) {
    abort_msg("Inconsistent window definition (begin or length) between components for data files");
  }
  if (strcmp(tmpstring,phase)) {
    sprintf(abort_str,"Inconsistent phase name between components for data files %s vs %s",tmpstring,phase);
   abort_msg(abort_str);
  } 
 
  if(par->method==CORREL) {
    ref1=find_window(par->method_q.cor_par.hdr_ref1,par->method_q.cor_par.data_ref1,
		     &par->window,0.0,&refwbeg,&refwlen,tmpstring);
    ref2=find_window(par->method_q.cor_par.hdr_ref2,par->method_q.cor_par.data_ref2,
		     &par->window,0.0,&ldum1,&ldum2,tmpstring);
    /* setting maxsumlag to 0.0 for reference trace as we only shift the data trace */
    if (ldum1 != refwbeg || ldum2 != refwlen) {
/*       fprintf(stderr,"DEBUG refwbeg %d ldum1 %d refwln %d ldum2 %d\n",refwbeg,ldum1,refwlen,ldum2); */
      abort_msg("Inconsistent window definition (begin or length) between components for reference files");
    }
    if ( refwlen != wlen ) {
      VRB(printf("refwlen %ld wlen %ld\n",refwlen,wlen));
      abort_msg("Length of analysis window for data and reference files must be identical");
    }
    if (strcmp(tmpstring,phase)) {
      abort_msg("Inconsistent phase names between data and reference files");
    } 
  }
  /* Initialise grid search */
  VRB(printf("Init grid search\n"));
  switch (par->model) {
  case SINGLE_HOR_SPLIT:
    hsplit=&(par->model_q.split_par.bot);
    m=(hsplit->fastmax-hsplit->fastmin+TOLERANCE)/hsplit->faststep + 1;
    n=(hsplit->timemax-hsplit->timemin+TOLERANCE)/hsplit->timestep + 1;
    VRB(printf("Grid search matrix size %d x %d\n",m,n));
    mat_res=gsl_matrix_alloc(m,n);
    if (par->method==MINEVALUE || par->method==MINTRANSVERSE ) 
      mat_pol=gsl_matrix_alloc(m,n);
    else if (par->method==CORREL) {
      mat_delay=gsl_matrix_alloc(m,n);
      mat_alpha=gsl_matrix_alloc(m,n);
    }
    break  ;
  case DOUBLE_HOR_SPLIT:
    abort_msg("Double layer splitting not implemented yet");
    break  ;
  }
  
  /* Single splitting implementation (later put this in loop over top layer for two layer splitting) */
  switch (par->method) {
  case MINEVALUE:
  case MINTRANSVERSE:
    mod_par=(par->method==MINEVALUE ? 3 : 2);
    single_split_sks(par->method, hsplit, hor1, hor2, wbeg, wlen, par->hdr_hor1->delta, par->hdr_hor1->baz, mat_res, mat_pol);
    break;
  case CONV:
    abort_msg("CONV method not implemented yet");
    break;
  case CORREL:
    mod_par=4;
    single_split_correl(par->method, hsplit, par->method_q.cor_par.maxshift, hor1, hor2, wbeg, ref1, ref2, refwbeg, wlen, par->hdr_hor1->delta, mat_res, mat_delay, mat_alpha);
    mat_pol=mat_delay; /* fiddle such that I can use standard err_single_split processing */
    break;
  }

  /* Now we have the error surface and have to pick the best value and do error analysis */
  if (par->model==SINGLE_HOR_SPLIT) {
    switch(par->method) {
    case MINEVALUE:
      /* Note: there are really three model parameters as the polarisation is undetermined 
       (none of the papers seem to do it this way, though, so maybe my thinking is wrong and
       it should be just 2 model parameters? */
      err_single_split_sks(par, mat_res, hor1,hor2,wbeg,wlen,1,phase, mat_pol, NULL, NULL, NULL, -1);
      break;
    case MINTRANSVERSE:
      err_single_split_sks(par, mat_res, hor1,hor2,wbeg,wlen,0,phase, NULL, NULL, NULL, NULL, -1);
      break;
    case CONV:
      abort_msg("CONV method not implemented yet");
      break;
    case CORREL:
      par->dof_s *=2;  /* multiply degrees of freedom per second by two as there are two traces! */
      err_single_split_sks(par, mat_res,hor1,hor2,wbeg,wlen,2,phase,mat_delay,mat_alpha,ref1,ref2,refwbeg);
      break;
    }
  }
  return(0);
}


void err_single_split_sks(ms_params *par, gsl_matrix *m_res_energy, gsl_vector_float *north, gsl_vector_float *east,long beg, long len, long mod_par, char *phase,  gsl_matrix *m_aux1, gsl_matrix *m_aux2, gsl_vector_float *ref_north, gsl_vector_float *ref_east, long refbeg) {
/* Error analysis and output 
   m_err     matrix containing residual energy estimates 
   len:      length of analysis window in samples
   mod_par:  number of additional model parameters (e.g. 0 for Min Transvers, or 1 for Min E Value (polarisation), 2 for reference station method (alpha value and delay)
   phase:    name of the phase to be used

Compact output format
Output:
summary of result to standard out

root_split.txt: one line summary
event station phase fast errfast(perp) split-delay errdelay(perp) reject-null? baz dist depth pol accept-null 

pol: inferred initial polarisation direction (either from BAZ or eigenvector)
accept-null  : 1 if pattern appears to be zero-splitting pattern 
               0 noisy or significant splitting (zero=0 and reject-null<0.95 implies noisy data)

root_err.hdr:   text file with information about _err.bin file
                  line0: file info

                  line1: dim      (number of dimensions of grid file)
		  line2: m n ...  (number of values along each dimension
		  line3: descr1 min1 max1 step1   data for first dimension
		  line4: descr2 min2 max2 step2   

		  Alternatives for line0
		  line0: MinTransverse mod_par dof
		  line0: MinEvalue     mod_par dof
		  line0: ResidualRefObsData     mod_par dof
		  line0: POLARIZATION  
		  (NB Only output polarisation file for MINEVALUE method)

root_err.bin:   Error surface as binary table

Optional: 
root_err.grd:   GMT format error surface (normalised such that minimum value = 1 but retain actual value in scale parameter
root_err.cont:   Contour file with 68%,95.2%,99%, 99.9%,99.99% and 99.999% confidence intervals 
*/
  double conf[9]={.68,  .95,.99, .999,.9999,.99999,.999999,.9999999, .99999999 };  /* conf[1] is the level of significance */
  double contour[9];
  double data_dof,dof;
  double null;
  double emin;
  double pol=NAN,lin=NAN;
  double best_time,best_fast;
  float delta=par->hdr_hor1->delta;
  float maxamp,tmpmax,tmpmin;
  double lbound_fast,ubound_fast,err_fast,lbound_time,ubound_time,err_time;
  size_t dum1,dum2;
  long split_par;
  long jpol=-1,jpolp=-1;
  long i,jmin,kmin,j1,j2,k1,k2,m,n,k1sec;
  hor_split *hsplit=&par->model_q.split_par.bot;
  int periodic,zero;
  char evname[256],tmpstring[256],dummystring[256],cmdstring[256],rejectstring[256];
  char *methodstring;
  char *descrip1=NULL,*ext1=NULL,*ext,*descrip;
  char *descrip2=NULL,*ext2=NULL;
  sachdr *hdr=par->hdr_hor1;
  FILE *output;
  int status;
  gsl_matrix *m_aux;
  gsl_vector_float_view vue1,vue2;
  gsl_vector_view dvue1;
  gsl_vector_float *vec1=gsl_vector_float_alloc(len);
  gsl_vector_float *vec2=gsl_vector_float_alloc(len);

  m=m_res_energy->size1;
  n=m_res_energy->size2;

  /* determine minimum solution */
  gsl_matrix_min_index(m_res_energy,&dum1,&dum2); jmin=(long)dum1; kmin=(long)dum2;
  emin=gsl_matrix_get(m_res_energy,jmin,kmin);
  best_fast=hsplit->fastmin+jmin*hsplit->faststep;
  best_time=hsplit->timemin+kmin*hsplit->timestep;

  split_par=2; /* Number of parameters of splitting model: Fast direction and splitting delay */
  /* Method specific values */
  switch (par->method) {
  case MINTRANSVERSE:
    methodstring=strdup("MinimumTransverse");
    pol=hdr->baz;
    /* polarity is BAZ mod 180 */
    if (pol<0) pol+=180;
    if (pol>180) pol-=180;
    break;
  case MINEVALUE:
    methodstring=strdup("MinimumEigenvalue");
    descrip1=strdup("Polarisation");
    ext1=strdup("pol");
    /* Get polarity from measurement */
    pol=gsl_matrix_get(m_aux1,jmin,kmin);
    break;
  case CORREL:
    methodstring=strdup("ResidualRefObsData");
    descrip1=strdup("Shift");
    ext1=strdup("dly");
    descrip2=strdup("AmplitudeFactor");
    ext2=strdup("alpha");

    vue1=gsl_vector_float_subvector(ref_north,refbeg,len);
    vue2=gsl_vector_float_subvector(ref_east,refbeg,len);
    pol=gsl_df_polarisation(&vue1.vector,&vue2.vector,&lin);
    break;
  }
  /* compute confidence levels (multipliers of minimum energy) */
  data_dof=par->dof_s * (len * delta - MAX(par->window.taper,0.0));
  dof=data_dof-mod_par-split_par;
  VRB(printf("data_dof %f mod_par %ld split_par %ld dof %f\n",data_dof,mod_par,split_par,dof));
  if (dof<2) {
    sprintf(warn_str,"Degrees of freedom less than 2 (DOF=%f). Set to 2 but error estimates likely to be meaningless.",dof);
    warn_msg(warn_str);
    dof=2;
  }
  if (!isnan(pol)) {
    /* index in table corresponding to inferred initial polarisation */
    jpol=ROUND((pol-hsplit->fastmin)/hsplit->faststep);  
    /* index in table corresponding to direction perpendicular to inferred initial polarisation */
    jpolp=ROUND(((pol>=90 ? pol-90 : pol+90)   -hsplit->fastmin)/hsplit->faststep);
  }
  /* index corresponding to splitting delay of 1s */
  k1sec=ROUND((1-hsplit->timemin)/hsplit->timestep);
  null=0.0;  

  /* check error bounds and likelyhood level of null splitting */
  for(i=0;i<9;i++) {
    contour[i]=1+split_par*invfisher((double)split_par,dof,conf[i])/dof;
    VRB(printf("i: %ld confidence %f contour %f\n",i,conf[i],contour[i]));
    if (F_EQ(hsplit->timemin,0.0) && contour[i]*emin<gsl_matrix_get(m_res_energy,0,0))
      null=conf[i];
  }
  if (par->make_grd) {
      strcpy(dummystring, par->root);
      strcat(dummystring, ".cont");
      output=fopen(dummystring,"w");
      if (!output) { 
	sprintf(abort_str,"Cannot open %s for output.",dummystring);
	abort_msg(abort_str);
      }
      for (i=0;i<9;i++){
	fprintf(output,"%f %s\n",contour[i]*emin, (i==1 ? "A" : "C" ));
      }
      fclose(output);
  }

  if(hsplit->timemin!=0.0)
    null=-1;    /* Null confidence level is only meaningful if first column represents 0 splitting */
  /* Test whether data is not only consistent with null splitting but  whether the patterns is the typical
     two-branch zero splitting pattern.  If
     1. fast or slow direction perpendicular to backazimuth and a splitting delay of 1s is accepted within 95% conf.
        interval
     2. fast or slow direction at 45 degree to backazimuth and a splitting delay of 1s is rejected at 95% confidence
        then call the result a likely null
	(NB the proper way to do this would be to compare the variances for these two splitting parameters with an 
	F-test)
  */
  VRB(printf("jpol: %ld jpolp %ld\n",jpol,jpolp));
  if ( jpol>=0 && jpol<m && jpolp>=0 && jpolp<m && k1sec>=0 && k1sec<n ) {
    if( contour[1]*emin>=gsl_matrix_get(m_res_energy,jpol,k1sec)
       && contour[1]*emin>=gsl_matrix_get(m_res_energy,jpolp,k1sec)
	&& contour[1]*emin<gsl_matrix_get(m_res_energy,(jpol+jpolp)/2,k1sec) )
      zero=1;    /* Zero splitting likely */
    else
      zero=0;    /* Zero splitting unlikely or can't tell */
  } else
    zero=-1;
  

  /* Splitting Time */
  /* search in perpendicular direction form minimum solution */
  for(k1=kmin+1;k1<n;k1++) {
    if (gsl_matrix_get(m_res_energy,jmin,k1)>contour[1]*emin) 
      break;
  }
  for(k2=kmin-1;k2>=0;k2--) {
    if (gsl_matrix_get(m_res_energy,jmin,k2)>contour[1]*emin)
      break;
  }
  lbound_time=hsplit->timemin+(k2+1)*hsplit->timestep;
  ubound_time=hsplit->timemin+(k1-1)*hsplit->timestep;
  err_time=MAX(ubound_time-best_time,best_time-lbound_time);
  VRB(printf("Local search time %f-%f (E: %f)\n",lbound_time,ubound_time,err_time));
  /* Global search, overwrite results of perpendicular search: */
  for (k1=n-1; k1>=kmin;k1--) {
    dvue1=gsl_matrix_column(m_res_energy,k1);
    if (gsl_vector_min(&dvue1.vector) <= contour[1]*emin) 
      break;
  }
  for (k2=0; k2<=kmin;k2++) {
    dvue1=gsl_matrix_column(m_res_energy,k2);
    if (gsl_vector_min(&dvue1.vector) <= contour[1]*emin) 
      break;
  }

  lbound_time=hsplit->timemin+k2*hsplit->timestep;
  ubound_time=hsplit->timemin+k1*hsplit->timestep;
  err_time=MAX(ubound_time-best_time,best_time-lbound_time);
  VRB(printf("Global search time %f-%f (E: %f)\n",lbound_time,ubound_time,err_time));
  if (k1==n-1 && k2==0)     /* both bounds at limit of grid search -> no constraints */
    /*     err_time=nan("");    */  /* set to NaN */ 
    err_time=NAN;
  else if (k1==n-1 )     /* upper bounds at limit of grid search */
    err_time=-err_time;    /* set error to negative (as flag) */

  /* Fast direction: there is a complication here because of the 180 deg periodicity of solutions.
     We resolve this by wrapping round the
     search but only if 180+fastmin-fastmax is less than 2 faststep (otherwise only a sub-region is searched
     and wrap around is not an issue */
  if(180+hsplit->fastmin-hsplit->fastmax<2*hsplit->faststep) {
    VRB(printf("Periodic fast direction\n"));
    periodic=1;
  } else {
    VRB(printf("Non-Periodic fast direction\n"));
    periodic=0; 
  }
  for(j1=jmin+1; j1<(periodic ? m/2+1+jmin: m); j1++) {
/*     VRB(printf("j1=%d Stop at %d\n",j1, (periodic ? jmin+m : m))); */
    if (gsl_matrix_get(m_res_energy,PER(m,j1),kmin)>contour[1]*emin) 
      break;
  }
/*   VRB(printf("Init j2=%d Stop at %d\n",jmin-1, (periodic ? jmin-m+1 : 0))); */
  for(j2=jmin-1; j2>=(periodic ? jmin-m/2 : 0) ;j2--) {
/*     VRB(printf("j2=%d Stop at %d\n",j2, (periodic ? jmin-m+1 : 0))); */
    if (gsl_matrix_get(m_res_energy,PER(m,j2),kmin)>contour[1]*emin) 
      break;
  }
  ubound_fast=hsplit->fastmin+(j1-1)*hsplit->faststep;
  lbound_fast=hsplit->fastmin+(j2+1)*hsplit->faststep;
  err_fast=MAX(ubound_fast-best_fast,best_fast-lbound_fast);
  VRB(printf("Local search fast direction %f-%f (E: %f)\n",lbound_fast,ubound_fast,err_fast));

  /* Global search, overwrite results of perpendicular search: */
  for (j1=(periodic ? jmin+m/2 : m-1); j1>=jmin;j1--) {
    dvue1=gsl_matrix_row(m_res_energy,PER(m,j1));
    if (gsl_vector_min(&dvue1.vector) <= contour[1]*emin) 
      break;
  }
  for (j2=(periodic ? jmin-m/2 : 0); j2<=jmin;j2++) {
    dvue1=gsl_matrix_row(m_res_energy,PER(m,j2));
    if (gsl_vector_min(&dvue1.vector) <= contour[1]*emin) 
      break;
  }
  ubound_fast=hsplit->fastmin+j1*hsplit->faststep;
  lbound_fast=hsplit->fastmin+j2*hsplit->faststep;
  err_fast=MAX(ubound_fast-best_fast,best_fast-lbound_fast);
  VRB(printf("Global search fast direction %f-%f (E: %f)\n",lbound_fast,ubound_fast,err_fast));

  if (err_fast>75.)
    err_fast=NAN;



  /* Output */

  /* stdout: */
  make_event_name(evname,hdr,EVN_YYJJJHHMM);
  if ((null)<0) 
    strcpy(rejectstring, "UNK");
  else if (null==0.0)
    strcpy(rejectstring, "<68");
  else
    sprintf(rejectstring,"%7.4f",null*100);
/*   strncpy(dummystring,hdr->kstnm,8); dummystring[8]='\0'; */
/*   VRB( printf("Energy:     %f\n",emin)); */

  
  printf("Event:            %s\n",evname); 
  printf("Station:          %8.8s\n",hdr->kstnm);
  printf("Phase:            %-8s\n",phase); 
  printf("%s: %f\n",methodstring,emin);
  printf("Fast direction:   %3.0f  +- %3.0f  (% 4.0f - %3.0f )\n",best_fast,err_fast,lbound_fast,ubound_fast); 
  printf("Splitting delay : %4.2f +- %4.2f (% 4.2f - %4.2f )\n",best_time,err_time,lbound_time,ubound_time); 
  if (strncmp(rejectstring,"UNK",3))
    printf("Reject Null (%%) : %s (E: %f)\n",rejectstring,gsl_matrix_get(m_res_energy,0,0));
  printf("Backazimuth (dg): %3.0f\n",hdr->baz);
  printf("Distance (dg):    %3.0f\n",hdr->gcarc);
  printf("Depth (km):       %3.0f\n",hdr->evdp);
  if (!isnan(pol))
    printf("Initial Pol.(dg): %3.0f\n",pol);
  if (!isnan(lin))
    printf("Linearity:        %4.2f\n",lin);
  if (descrip1) 
    printf("%-18s%f\n",descrip1,gsl_matrix_get(m_aux1,jmin,kmin));
  if (descrip2) 
    printf("%-18s%f\n",descrip2,gsl_matrix_get(m_aux2,jmin,kmin));
  VRB(printf("Zero: %d, null %f\n",zero,null));
  VRB(printf("Null %f conf[1] %f cond %d\n",null,conf[1],null>=conf[1]));

  if(zero==1) {
    if (null<conf[1]){
      printf("Zero splitting or aligned fast direction likely\n");
    }  else {
      printf("Null pattern but null-splitting rejected: Data likely to be noisy or splitting model not sufficiently complex\n");
    }
  } else if (zero==0) {
    if (null>=conf[1]){
      printf("Significant splitting likely\n");
    }  else {
      printf("Null accepted but no null pattern: Data likely to be noisy or splitting model not sufficiently complex\n");
    }
  }

  /* compact one line output */
  strcpy(dummystring, par->root);
  strcat(dummystring, "_split.txt");
  output=fopen(dummystring,"w");
  if (!output) { 
    sprintf(abort_str,"Cannot open %s for output.",dummystring);
    abort_msg(abort_str);
  }

  fprintf(output,
	  "#Event       Station  Phase   Fast Err SDel Err   Rej-0   BAZ Dis Dpt Pol Emin   %-7s %-7s Zero?\n",
	  ext1 ? ext1 : "", ext2 ? ext2 : "");

  fprintf(output,
	  "%s %8.8s %-8s %3.0f %3.0f %4.2f %5.2f %-7.7s %3.0f %3.0f %3.0f %3.0f %6.4f %-7.3f %-7.3f %1s\n",
	  evname, 
	  hdr->kstnm,
	  phase,
	  best_fast,err_fast,
	  best_time,err_time,
	  rejectstring,
	  hdr->baz,
	  hdr->gcarc,
	  hdr->evdp,
	  pol,
	  emin,
	  m_aux1 ? gsl_matrix_get(m_aux1,jmin,kmin) : NAN,
	  m_aux2 ? gsl_matrix_get(m_aux2,jmin,kmin) : NAN,
	  zero==-1 ? "?" : (zero==1?"1":"0") );
  fclose(output);

  /* Err_hdr */
  strcpy(tmpstring, par->root);
  strcat(tmpstring, "_err.hdr");
  output=fopen(tmpstring,"w");
  if (!output) { 
    sprintf(abort_str,"Cannot open %s for output.",tmpstring);
    abort_msg(abort_str);
  }
  fprintf(output,"%s %ld %f\n",methodstring,split_par,dof);
  fprintf(output,"2\n");
  fprintf(output,"%ld %ld\n",m,n);
  fprintf(output,"Fast     %f %f %f\n",hsplit->fastmin,hsplit->fastmin+(m-1)*hsplit->faststep,hsplit->faststep);
  fprintf(output,"SplitDly %f %f %f\n",hsplit->timemin,hsplit->timemin+(n-1)*hsplit->timestep,hsplit->timestep);
  fclose(output);
  /* err bin */
  strcpy(tmpstring, par->root);
  strcat(tmpstring, "_err.bin");
  output=fopen(tmpstring,"wb");
  if (!output) { 
    sprintf(abort_str,"Cannot open %s for output.",tmpstring);
    abort_msg(abort_str);
  }
  gsl_matrix_fwrite(output, m_res_energy);
  fclose(output);

  if (par->make_grd ) {
    sprintf(cmdstring,"xyz2grd %s_err.bin -DSplitDly/Fast/%s/1/0/SingleSplit/\"Created by multisplit\" -G%s_err.grd -I%f/%f -R%f/%f/%f/%f -ZBLd",
	    par->root,methodstring,par->root,
	    hsplit->timestep,hsplit->faststep,
	    hsplit->timemin,hsplit->timemin+(n-1)*hsplit->timestep,hsplit->fastmin,hsplit->fastmin+(m-1)*hsplit->faststep);
    VRB(printf("GRD conversion: %s\n",cmdstring));
    status=system(cmdstring);
    if (status) {
      sprintf(warn_str,"External GMT command xyz2grd execution failed. Error status: %d",status);
      warn_msg(warn_str);
    }
  }

  /* treat both auxilary variables the same way */
  for(i=0; i<2;i++) {
    switch (i) {
    case 0:
      m_aux=m_aux1;
      descrip=descrip1;
      ext=ext1;
      break;
    case 1:
      m_aux=m_aux2;
      descrip=descrip2;
      ext=ext2;
      break;
    }
    VRB(printf("i %ld M_AUX %p %p %p\n",i,m_aux,m_aux1,m_aux2));

    if ( m_aux && descrip ) {
      VRB(printf("Generating bin files for %s\n",descrip));
      sprintf(tmpstring,"%s_%s.hdr",par->root,ext);
      /*     strcpy(tmpstring, par->root); */
      /*     strcat(tmpstring, "_pol.hdr"); */
      output=fopen(tmpstring,"w");
      if (!output) { 
	sprintf(abort_str,"Cannot open %s for output.",tmpstring);
	abort_msg(abort_str);
      }
      fprintf(output,"%s\n",descrip);
      fprintf(output,"2\n");
      fprintf(output,"%ld %ld\n",m,n);
      fprintf(output,"Fast     %f %f %f\n",hsplit->fastmin,hsplit->fastmin+(m-1)*hsplit->faststep,hsplit->faststep);
      fprintf(output,"SplitDly %f %f %f\n",hsplit->timemin,hsplit->timemin+(n-1)*hsplit->timestep,hsplit->timestep);
      fclose(output);
      /* err bin */
      sprintf(tmpstring,"%s_%s.bin",par->root,ext);
      /*     strcpy(tmpstring, par->root); */
      /*     strcat(tmpstring, "_pol.bin"); */
      output=fopen(tmpstring,"wb");
      if (!output) { 
	sprintf(abort_str,"Cannot open %s for output.",tmpstring);
	abort_msg(abort_str);
      }
      gsl_matrix_fwrite(output, m_aux);
      fclose(output);
      if (par->make_grd ) {
	sprintf(cmdstring,"xyz2grd %s_%s.bin -DSplitDly/Fast/%s/1/0/SingleSplit/\"Created by multisplit - minevalue mode\" -G%s_%s.grd -I%f/%f -R%f/%f/%f/%f -ZBLd",
		par->root,ext,descrip,par->root,ext,
		hsplit->timestep,hsplit->faststep,
		hsplit->timemin,hsplit->timemin+(n-1)*hsplit->timestep,hsplit->fastmin,hsplit->fastmin+(m-1)*hsplit->faststep);
	VRB(printf("GRD conversion: %s\n",cmdstring));
	status=system(cmdstring);
	if (status) {
	  sprintf(warn_str,"External GMT command xyz2grd execution failed. Error status: %d",status);
	  warn_msg(warn_str);
	}
      }
    }
  }
  if ( par->make_grd | MAKE_GMT  ) {
    /* Output time series */
    long itime;

    vue1=gsl_vector_float_subvector(north,beg,len);
    vue2=gsl_vector_float_subvector(east,beg,len);
    
    gsl_vector_float_memcpy(vec1, &vue1.vector);
    gsl_vector_float_memcpy(vec2, &vue2.vector);

    maxamp=0;
/*    output=open_for_write(par->root,"_north.xy"); */
/*     gsl_float_write_timeseries(output,vec1,0,delta); */
/*     fclose(output); */

/*    output=open_for_write(par->root,"_east.xy"); */
/*     gsl_float_write_timeseries(output,vec2,0,delta); */
/*     fclose(output); */

    /* Rotate to baz */
    gsl_float_rotate(vec1,vec2,180-hdr->baz);
 
    VRB(printf("Writing radial\n"));
    gsl_vector_float_minmax(vec1,&tmpmin,&tmpmax);
    maxamp=MAX(maxamp,MAX(tmpmax,-tmpmin));
    output=open_for_write(par->root,"_rad.xy");
    gsl_float_write_timeseries(output,vec1,0,delta);
    fclose(output);

    VRB(printf("Writing transverse\n"));
    gsl_vector_float_minmax(vec2,&tmpmin,&tmpmax);
    maxamp=MAX(maxamp,MAX(tmpmax,-tmpmin));
    output=open_for_write(par->root,"_tra.xy");
    gsl_float_write_timeseries(output,vec2,0,delta);
    fclose(output);


/*     VRB(printf("Writing rad/tra particle motion\n")); */
/*     output=open_for_write(par->root,"_rt_pmp.xy"); */
/*     gsl_float_write_pmp(output,vec2,vec1); */
/*     fclose(output); */

    /* rotate to F/S */
    gsl_float_rotate(vec1,vec2,180+hdr->baz-best_fast);
    
/*     gsl_vector_float_memcpy(vec1, north); */
/*     gsl_vector_float_memcpy(vec2, east); */
/*     gsl_float_rotate(vec1,vec2,-hdr->baz+180); */
 
    VRB(printf("Writing fast\n"));
    gsl_vector_float_minmax(vec1,&tmpmin,&tmpmax);
    maxamp=MAX(maxamp,MAX(tmpmax,-tmpmin));
    output=open_for_write(par->root,"_fast.xy");
    gsl_float_write_timeseries(output,vec1,0,delta);
    fclose(output);

    VRB(printf("Writing slow\n"));
    gsl_vector_float_minmax(vec2,&tmpmin,&tmpmax);
    maxamp=MAX(maxamp,MAX(tmpmax,-tmpmin));
    output=open_for_write(par->root,"_slow.xy");
    gsl_float_write_timeseries(output,vec2,0,delta);
    fclose(output);

/*     VRB(printf("Writing fast/slow particle motion\n")); */
/*     output=open_for_write(par->root,"_fs_pmp.xy"); */
/*     gsl_float_write_pmp(output,vec2,vec1); */
/*     fclose(output); */

    itime=(long)ROUND(best_time/delta);
    vue1=gsl_vector_float_subvector(vec1,0,len-itime); /* retarded fast */
    vue2=gsl_vector_float_subvector(vec2,itime,len-itime); /* time advanced slow */

    VRB(printf("Writing fastcor\n"));
    gsl_vector_float_minmax(&vue1.vector,&tmpmin,&tmpmax);
    maxamp=MAX(maxamp,MAX(tmpmax,-tmpmin));
    output=open_for_write(par->root,"_fastcor.xy");
    gsl_float_write_timeseries(output,&(vue1.vector),(itime/2)*delta,delta);
    fclose(output);

    VRB(printf("Writing slowcor\n"));
    gsl_vector_float_minmax(&vue2.vector,&tmpmin,&tmpmax);
    maxamp=MAX(maxamp,MAX(tmpmax,-tmpmin));
    output=open_for_write(par->root,"_slowcor.xy");
    gsl_float_write_timeseries(output,&(vue2.vector),(itime/2)*delta,delta);
    fclose(output);

/*     VRB(printf("Writing fast/slowcor particle motion\n")); */
/*     output=open_for_write(par->root,"_fscor_pmp.xy"); */
/*     gsl_float_write_pmp(output,&vue2.vector,&vue1.vector); */
/*     fclose(output); */

    /* rotate back to BAZ */
    gsl_float_rotate(&vue1.vector,&vue2.vector,180-hdr->baz+best_fast);

    VRB(printf("Writing radcor\n"));
    gsl_vector_float_minmax(&vue1.vector,&tmpmin,&tmpmax);
    maxamp=MAX(maxamp,MAX(tmpmax,-tmpmin));
    output=open_for_write(par->root,"_radcor.xy");
    gsl_float_write_timeseries(output,&(vue1.vector),(itime/2)*delta,delta);
    fclose(output);

    VRB(printf("Writing tracor\n"));
    gsl_vector_float_minmax(&vue2.vector,&tmpmin,&tmpmax);
    maxamp=MAX(maxamp,MAX(tmpmax,-tmpmin));
    output=open_for_write(par->root,"_tracor.xy");
    gsl_float_write_timeseries(output,&(vue2.vector),(itime/2)*delta,delta);
    fclose(output);

/*     VRB(printf("Writing rad/tracor particle motion\n")); */
/*     output=open_for_write(par->root,"_rtcor_pmp.xy"); */
/*     gsl_float_write_pmp(output,&vue2.vector,&vue1.vector); */
/*     fclose(output); */

    VRB(printf("Writing GMT script\n"));
    output=open_for_write(par->root,".gmt");
    fprintf(output,"#/bin/csh\n");
    fprintf(output,"# script auto-generated by multisplit\n");
    fprintf(output,"\n");
    fprintf(output,"# Station-event dependent:\n");
    fprintf(output,"set root=%s\n",par->root);
    fprintf(output,"set maxamp=%g\n",maxamp);
    fprintf(output,"set bestfast=%f\n",best_fast);
    fprintf(output,"set besttime=%f\n",best_time);
    fprintf(output,"set baz=%f\n",hdr->baz);
    if (par->make_grd == MAKE_GMT5) {
      fprintf(output,"alias grdinfo gmt grdinfo\n");
      fprintf(output,"alias psxy gmt psxy\n");
      /* continue here to put gmt5 commands */
    }
    if(!isnan(pol)) 
      fprintf(output,"set pol=%f\n",pol);
    fprintf(output,"cat > $root.description <<EOF\n");
    fprintf(output,"> 10 29 14 0 0 CT 0.564 20  c\n");
    fprintf(output,"%s %8.8s %-8s BAZ %3.0f Dist %3.0f Dp %3.0f\n\n",evname,hdr->kstnm,phase,hdr->baz, hdr->gcarc, hdr->evdp);
    fprintf(output,"%s: Res %f",methodstring,emin);
    if (descrip1 && m_aux1)
      fprintf(output," %s %f",descrip1,gsl_matrix_get(m_aux1,jmin,kmin));
    if (descrip2 && m_aux2)
      fprintf(output," %s %f",descrip2,gsl_matrix_get(m_aux2,jmin,kmin));
    fprintf(output," \n\nFast %3.0f \\261 %3.0f Time %4.2f \\261 %4.2f RejectNull %s",
	    best_fast,err_fast, best_time,err_time,rejectstring);
    if (!isnan(pol))  
      fprintf(output," Pol %3.0f (Dev: %ld)\n",pol,PER(180,(long)(pol-hdr->baz+0.5+90))-90 );
    else 
      fprintf(output,"\n");
    fprintf(output,"EOF\n\n");
    fprintf(output,"\
### Everything below this line is independent of the particular event used\n\
set grdrange=`grdinfo -C ${root}_err.grd | awk '{print $2 \"/\" $3 \"/\" $4 \"/\" $5 }'`\n\
set maxamp2=`echo $maxamp | awk '{ print 2*$1}'`\n\
\n\
set timerange=( `awk 'NR==1 { print $1 } { lastx=$1 } END { print lastx }' ${root}_rad.xy` )\n\
\n\
set psfile=${root}.ps\n\
\n\
gmtdefaults -D >.gmtdefaults\n\
gmtset PAGE_ORIENTATION portrait MEASURE_UNIT cm WANT_EURO_FONT TRUE LABEL_FONT_SIZE 12 ANOT_FONT_SIZE 10 PAPER_MEDIA a4 D_FORMAT %%lg\n\
# 3cm Descriptive text\n\
pstext -M -X0 -Y0 -R0/20/0/29 -Jx1 -K > $psfile <${root}.description\n\
\n\
# 8 cm Error surface\n\
grdcontour -X2 -Y20.5 ${root}_err.grd -C${root}.cont -R$grdrange -JX17/6.5 -B0.5:\"Splitting Delay (s)\":/20:\"Fast direction\":WSen -O -K -A-1f1 -G1000 -Wa1.5p -Wc0.5p >>$psfile\n\
psxy -R -JX -Sx0.3 -W2p  -O -K  >>$psfile <<EOF\n\
$besttime $bestfast\n\
EOF\n\
psxy -R -JX -W1p/200/200/200to -O -K  >>$psfile <<EOF\n\
0 $baz\n\
10 $baz\n\
EOF\n\
if ( $?pol ) then \n\
  set polp=`echo $pol | awk '{ if ($1<90) {print $1+90;} else print $1-90 }'`\n\
  psxy -R -JX -M -W1p/200/200/200 -O -K  >>$psfile <<EOF\n\
>\n\
0 $pol\n\
10 $pol\n\
> \n\
0 $polp\n\
10 $polp\n\
EOF\n\
endif\n\
# 4 cm Rad Transverse\n\
psxy ${root}_rad.xy -Y-6.5 -R$timerange[1]/$timerange[2]/-$maxamp/$maxamp -JX13/4 -Ba5f1/${maxamp}::wseN -W1p -O -K >>$psfile\n\
psxy ${root}_tra.xy -R$timerange[1]/$timerange[2]/-$maxamp/$maxamp -JX -W1p0/0/0to -O -K >>$psfile\n\
pstext <<EOF -JX -R0/1/0/1 -O -K >>$psfile\n\
0.05 0.95 12 0 0 LT Radial-Transverse\n\
EOF\n\
paste ${root}_rad.xy ${root}_tra.xy | awk '{ print $4,$2 }' | psxy -X13 -R-${maxamp}/${maxamp}/-${maxamp}/${maxamp} -JX4/4 -W1p -B${maxamp}:\"Transverse\":/${maxamp2}:\"Radial\":wsen  -O -K >>$psfile\n\
#psxy ${root}_rt_pmp.xy -X13 -R-${maxamp}/${maxamp}/-${maxamp}/${maxamp} -JX4/4 -W1p -B${maxamp}:\"Transverse\":/${maxamp2}:\"Radial\":wsen  -O -K >>$psfile\n\
\n\
# 4 cm Fast Slow\n\
psxy ${root}_fast.xy -X-13 -Y-4 -R$timerange[1]/$timerange[2]/-$maxamp/$maxamp -JX13/4 -Ba5f1/${maxamp}::wsen -W1p -O -K >>$psfile\n\
psxy ${root}_slow.xy -R$timerange[1]/$timerange[2]/-$maxamp/$maxamp -JX -W1p0/0/0to -O -K >>$psfile\n\
pstext <<EOF -JX -R0/1/0/1 -O -K >>$psfile\n\
0.05 0.95 12 0 0 LT Fast-Slow\n\
EOF\n\
paste ${root}_fast.xy ${root}_slow.xy | awk '{ print $4,$2 }' | psxy -X13  -R-${maxamp}/${maxamp}/-${maxamp}/${maxamp} -JX4/4 -W1p -B${maxamp}:\"Slow\":/${maxamp2}:\"Fast\":wsen  -O -K >>$psfile\n\
#psxy ${root}_fs_pmp.xy -X13  -R-${maxamp}/${maxamp}/-${maxamp}/${maxamp} -JX4/4 -W1p -B${maxamp}:\"Slow\":/${maxamp2}:\"Fast\":wsen  -O -K >>$psfile\n\
\n\
# 4 cm Fast Slow, corrected\n\
psxy ${root}_fastcor.xy -X-13 -Y-4 -R$timerange[1]/$timerange[2]/-$maxamp/$maxamp -JX13/4 -Ba5f1/${maxamp}::wsen -W1p -O -K >>$psfile\n\
psxy ${root}_slowcor.xy -R$timerange[1]/$timerange[2]/-$maxamp/$maxamp -JX -W1p0/0/0to -O -K >>$psfile\n\
pstext <<EOF -JX -R0/1/0/1 -O -K >>$psfile\n\
0.05 0.95 12 0 0 LT Fast-Slow, corrected\n\
EOF\n\
paste ${root}_fastcor.xy ${root}_slowcor.xy | awk '{ print $4,$2 }' | psxy -X13  -R-${maxamp}/${maxamp}/-${maxamp}/${maxamp} -JX4/4 -W1p -B${maxamp}:\"Slow\":/${maxamp2}:\"Fast\":wsen  -O -K >>$psfile\n\
#psxy ${root}_fscor_pmp.xy -X13  -R-${maxamp}/${maxamp}/-${maxamp}/${maxamp} -JX4/4 -W1p -B${maxamp}:\"Slow\":/${maxamp2}:\"Fast\":wsen  -O -K >>$psfile\n\
\n\
\n\
# 4 cm Rad Transverse, corrected \n\
psxy ${root}_radcor.xy -X-13 -Y-4 -R$timerange[1]/$timerange[2]/-$maxamp/$maxamp -JX13/4 -Ba5f1:\"Time (s)\":/${maxamp}::wSen -W1p -O -K >>$psfile\n\
psxy ${root}_tracor.xy -R$timerange[1]/$timerange[2]/-$maxamp/$maxamp -JX -W1p0/0/0to -O -K >>$psfile\n\
pstext <<EOF -JX -R0/1/0/1 -O -K >>$psfile\n\
0.05 0.95 12 0 0 LT Radial-Transverse, corrected\n\
EOF\n\
paste ${root}_radcor.xy ${root}_tracor.xy | awk '{ print $4,$2 }' | psxy -X13 -R-${maxamp}/${maxamp}/-${maxamp}/${maxamp} -JX4/4 -W1p -B${maxamp}:\"Transverse\":/${maxamp2}:\"Radial\":wsen  -O -K >>$psfile\n\
#psxy ${root}_rtcor_pmp.xy -X13 -R-${maxamp}/${maxamp}/-${maxamp}/${maxamp} -JX4/4 -W1p -B${maxamp}:\"Transverse\":/${maxamp2}:\"Radial\":wsen  -O -K >>$psfile\n\
\n\
\n\
psxy < /dev/null -Jx1 -R -O >>$psfile\n\
");
    fclose(output);

    VRB(printf("Executing GMT script\n"));
    sprintf(cmdstring,"csh %s.gmt",par->root);
    status=system(cmdstring);
    if (status) {
      sprintf(warn_str,"Execution of GMT script failed. Error status: %d",status);
      warn_msg(warn_str);
    }
    if ( ref_north && ref_east && !strcmp(ext1,"dly")) {
      double delay=gsl_matrix_get(m_aux1,jmin,kmin);
      double alpha=gsl_matrix_get(m_aux2,jmin,kmin);

      vue1=gsl_vector_float_subvector(ref_north,refbeg,len);
      vue2=gsl_vector_float_subvector(ref_east,refbeg,len);
      gsl_vector_float_memcpy(vec1, &vue1.vector);
      gsl_vector_float_memcpy(vec2, &vue2.vector);
      gsl_float_rotate(vec1,vec2,180-hdr->baz);

      VRB(printf("Writing radial reference trace\n"));
      gsl_vector_float_minmax(vec1,&tmpmin,&tmpmax);
      maxamp=MAX(maxamp,MAX(tmpmax,-tmpmin));
      output=open_for_write(par->root,"_refrad.xy");
      gsl_float_write_timeseries(output,vec1,0,delta);
      fclose(output);

      VRB(printf("Writing transverse reference trace\n"));
      gsl_vector_float_minmax(vec2,&tmpmin,&tmpmax);
      maxamp=MAX(maxamp,MAX(tmpmax,-tmpmin));
      output=open_for_write(par->root,"_reftra.xy");
      gsl_float_write_timeseries(output,vec2,0,delta);
      fclose(output);

      output=open_for_write(par->root,"-aux.gmt");
      fprintf(output,"#/bin/csh\n");
      fprintf(output,"# script auto-generated by multisplit\n");
      fprintf(output,"\n");
      fprintf(output,"# Station-event dependent:\n");
      fprintf(output,"set root=%s\n",par->root);
      fprintf(output,"set maxamp=%g\n",maxamp);
      fprintf(output,"set bestfast=%f\n",best_fast);
      fprintf(output,"set besttime=%f\n",best_time);
      fprintf(output,"set delay=%f\n",delay);
      fprintf(output,"set alpha=%f\n",alpha);
      fprintf(output,"set baz=%f\n",hdr->baz);
      fprintf(output,"cat > $root.description <<EOF\n");
      fprintf(output,"> 10 29 14 0 0 CT 0.564 20  c\n");
      fprintf(output,"%s %8.8s %-8s BAZ %3.0f Dist %3.0f Dp %3.0f\n\n",evname,hdr->kstnm,phase,hdr->baz, hdr->gcarc, hdr->evdp);
      fprintf(output,"%s: Res %f",methodstring,emin);
      if (descrip1 && m_aux1)
	fprintf(output,"%s %f ",descrip1,gsl_matrix_get(m_aux1,jmin,kmin));
      if (descrip2 && m_aux2)
	fprintf(output,"%s %f",descrip2,gsl_matrix_get(m_aux2,jmin,kmin));
      fprintf(output," \n\nFast %3.0f \\261 %3.0f Time %4.2f \\261 %4.2f RejectNull %s",
	      best_fast,err_fast, best_time,err_time,rejectstring);
      if (!isnan(pol))  
	fprintf(output," Pol %3.0f (Dev: %ld)\n",pol,PER(180,(long)(pol-hdr->baz+0.5+90))-90 );
      else 
	fprintf(output,"\n");
      fprintf(output,"EOF\n\n"); 

      fprintf(output,"\
### Everything below this line is independent of the particular event used\n\
\n\
set grdrange=`grdinfo -C ${root}_dly.grd | awk '{print $2 \"/\" $3 \"/\" $4 \"/\" $5 }'`\n\
set maxamp2=`echo $maxamp | awk '{ print 2*$1}'`\n\
\n\
set timerange=( `awk 'NR==1 { print $1 } { lastx=$1 } END { print lastx }' ${root}_rad.xy` )\n\
\n\
set psfile=${root}-aux.ps\n\
\n\
gmtdefaults -D >.gmtdefaults\n\
gmtset PAGE_ORIENTATION portrait MEASURE_UNIT cm WANT_EURO_FONT TRUE LABEL_FONT_SIZE 12 ANOT_FONT_SIZE 10 PAPER_MEDIA a4 D_FORMAT %%lg\n\
# 3cm Descriptive text\n\
pstext -M -X0 -Y0 -R0/20/0/29 -Jx1 -K > $psfile <${root}.description\n\
\n\
# Map delay\n\
grdcontour -X2 -Y20.5 ${root}_dly.grd -C0.1 -A0.5f7 -R$grdrange -JX17/6.5 -B0.5:\"Splitting Delay (s)\":/20:\"Fast direction\":WSen -O -K -G4c -Wa1.5p -Wc0.5p >>$psfile\n\
psxy -R -JX -Sx0.3 -W2p  -O -K  >>$psfile <<EOF\n\
$besttime $bestfast\n\
EOF\n\
# overlay 95%% confidence contour in gray\n\
awk '$2 == \"A\" { print  }' ${root}.cont > ${root}_tmp.acont\n\
grdcontour ${root}_err.grd -C${root}_tmp.acont -M -D${root}_tmp.95cont -JX  -R > /dev/null\n\
psxy ${root}_tmp.95cont -M  -R -JX -W1p/200/200/200ta  -O -K  >>$psfile \n\
\n\
# 4 cm RefRad refTransverse\n\
psxy ${root}_refrad.xy -Y-6.5 -R$timerange[1]/$timerange[2]/-$maxamp/$maxamp -JX13/4 -Ba5f1/${maxamp}::wseN -W1p -O -K >>$psfile\n\
psxy ${root}_reftra.xy -R$timerange[1]/$timerange[2]/-$maxamp/$maxamp -JX -W1p0/0/0to -O -K >>$psfile\n\
pstext <<EOF -JX -R0/1/0/1 -O -K >>$psfile\n\
0.05 0.95 12 0 0 LT Reference Radial-Transverse\n\
EOF\n\
paste ${root}_refrad.xy ${root}_reftra.xy | awk '{ print $4,$2 }' | psxy -X13 -R-${maxamp}/${maxamp}/-${maxamp}/${maxamp} -JX4/4 -W1p -B${maxamp}:\"Transverse\":/${maxamp2}:\"Radial\":wsen  -O -K >>$psfile\n\
\n\
# 4 cm Data radial-transverse\n\
psxy ${root}_rad.xy -X-13 -Y-4 -R$timerange[1]/$timerange[2]/-$maxamp/$maxamp -JX13/4 -Ba5f1/${maxamp}::wsen -W1p -O -K >>$psfile\n\
psxy ${root}_tra.xy -R$timerange[1]/$timerange[2]/-$maxamp/$maxamp -JX -W1p0/0/0to -O -K >>$psfile\n\
pstext <<EOF -JX -R0/1/0/1 -O -K >>$psfile\n\
0.05 0.95 12 0 0 LT Observed Radial-Transverse\n\
EOF\n\
paste ${root}_rad.xy ${root}_tra.xy | awk '{ print $4,$2 }' | psxy -X13  -R-${maxamp}/${maxamp}/-${maxamp}/${maxamp} -JX4/4 -W1p -B${maxamp}:\"Transverse\":/${maxamp2}:\"Radial\":wsen  -O -K >>$psfile\n\
\n\
# make amplitude and delay corrected traces\n\
awk '{ print $1-delay,$2*alpha/(1-alpha) }' delay=$delay alpha=$alpha ${root}_radcor.xy > ${root}_tmp_radcorcor.xy\n\
awk '{ print $1-delay,$2*alpha/(1-alpha) }' delay=$delay alpha=$alpha ${root}_tracor.xy > ${root}_tmp_tracorcor.xy\n\
# Interpolate reference station to grid of those data\n\
sample1d  ${root}_refrad.xy -N${root}_tmp_radcorcor.xy > ${root}_tmp_refrad.xy\n\
sample1d  ${root}_reftra.xy -N${root}_tmp_tracorcor.xy > ${root}_tmp_reftra.xy\n\
# Resample corrected data again on this grid to make sure we really share the same grid\n\
sample1d  ${root}_tmp_radcorcor.xy -N${root}_tmp_refrad.xy > ${root}_tmp_radcc.xy\n\
sample1d  ${root}_tmp_tracorcor.xy -N${root}_tmp_reftra.xy > ${root}_tmp_tracc.xy\n\
# Residual Observed-reference\n\
paste ${root}_tmp_radcc.xy ${root}_tmp_refrad.xy | awk '{ print $1,$2-$4 }' > ${root}_tmp_resrad.xy\n\
paste ${root}_tmp_tracc.xy ${root}_tmp_reftra.xy | awk '{ print $1,$2-$4 }' > ${root}_tmp_restra.xy\n\
\n\
# 4 cm  radial transverse, corrected (including amplitude and delay correction)\n\
psxy ${root}_tmp_radcorcor.xy -X-13 -Y-4 -R$timerange[1]/$timerange[2]/-$maxamp/$maxamp -JX13/4 -Ba5f1/${maxamp}::wsen -W1p -O -K >>$psfile\n\
psxy ${root}_tmp_tracorcor.xy -R$timerange[1]/$timerange[2]/-$maxamp/$maxamp -JX -W1p0/0/0to -O -K >>$psfile\n\
pstext <<EOF -JX -R0/1/0/1 -O -K >>$psfile\n\
0.05 0.95 12 0 0 LT Corrected Radial-Transverse\n\
EOF\n\
paste ${root}_tmp_radcorcor.xy ${root}_tmp_tracorcor.xy  | awk '{ print $4,$2 }' | psxy -X13  -R-${maxamp}/${maxamp}/-${maxamp}/${maxamp} -JX4/4 -W1p -B${maxamp}:\"Transverse\":/${maxamp2}:\"Radial\":wsen  -O -K >>$psfile\n\
\n\
\n\
# 4 cm Residual (Observed - reference)\n\
psxy ${root}_tmp_resrad.xy -X-13 -Y-4 -R$timerange[1]/$timerange[2]/-$maxamp/$maxamp -JX13/4 -Ba5f1:\"Time (s)\":/${maxamp}::wSen -W1p -O -K >>$psfile\n\
psxy ${root}_tmp_restra.xy -R$timerange[1]/$timerange[2]/-$maxamp/$maxamp -JX -W1p0/0/0to -O -K >>$psfile\n\
pstext <<EOF -JX -R0/1/0/1 -O -K >>$psfile\n\
0.05 0.95 12 0 0 LT Residual Radial-Transverse\n\
EOF\n\
paste ${root}_tmp_resrad.xy ${root}_tmp_restra.xy | awk '{ print $4,$2 }' | psxy -X13 -R-${maxamp}/${maxamp}/-${maxamp}/${maxamp} -JX4/4 -W1p -B${maxamp}:\"Transverse\":/${maxamp2}:\"Radial\":wsen  -O -K >>$psfile\n\
#psxy ${root}_rtcor_pmp.xy -X13 -R-${maxamp}/${maxamp}/-${maxamp}/${maxamp} -JX4/4 -W1p -B${maxamp}:\"Transverse\":/${maxamp2}:\"Radial\":wsen  -O -K >>$psfile\n\
\n\
psxy < /dev/null -Jx1 -R -O >>$psfile\n\
\\rm ${root}_tmp*\n\
"); 
      fclose(output);     
      VRB(printf("Executing GMT-aux script\n")); 
      sprintf(cmdstring,"csh %s-aux.gmt",par->root);
      status=system(cmdstring);
      if (status) {
	sprintf(warn_str,"Execution of GMT script failed. Error status: %d",status);
	warn_msg(warn_str);
      }
      VRB(printf("Writing GMT script colour\n"));
      output=open_for_write(par->root,"-colour.gmt");
      fprintf(output,"#/bin/csh\n");
      fprintf(output,"# script auto-generated by multisplit\n");
      fprintf(output,"\n");
      fprintf(output,"# Station-event dependent:\n");
      fprintf(output,"set root=%s\n",par->root);
      fprintf(output,"set maxamp=%g\n",maxamp);
      fprintf(output,"set bestfast=%f\n",best_fast);
      fprintf(output,"set besttime=%f\n",best_time);
      fprintf(output,"set delay=%f\n",delay);
      fprintf(output,"set alpha=%f\n",alpha);
      fprintf(output,"set baz=%f\n",hdr->baz);
      fprintf(output,"cat > $root.description <<EOF\n");
      fprintf(output,"> 5.7 29.5 12 0 0 CT 0.564 20  c\n");
      fprintf(output,"%s %8.8s %-8s BAZ %3.0f\\232 Dist %3.0f\\232 Dep %3.0f\n",evname,hdr->kstnm,phase,hdr->baz, hdr->gcarc, hdr->evdp);
      /*     fprintf(output,"%s: Res %f",methodstring,emin); */
      fprintf(output,"\nFast %3.0f \\261 %3.0f SplittingDelay %4.2f \\261 %4.2f\n\n",
	      best_fast,err_fast, best_time,fabs(err_time));
      if (descrip1 && m_aux1)
	fprintf(output,"%s %.2f ",descrip1,gsl_matrix_get(m_aux1,jmin,kmin));
      if (descrip2 && m_aux2)
	fprintf(output,"%s %.2f",descrip2,gsl_matrix_get(m_aux2,jmin,kmin));
      if (!isnan(pol) && !m_aux1)  
	fprintf(output," Pol %3.0f (Dev: %ld)\n",pol,PER(180,(long)(pol-hdr->baz+0.5+90))-90 );
      else 
	fprintf(output,"\n");
      fprintf(output,"EOF\n\n");
      fprintf(output,"\
### Everything below this line is independent of the particular event used\n\
set grdrange=`grdinfo -C ${root}_err.grd | awk '{print $2 \"/\" $3 \"/\" $4 \"/\" $5 }'`\n\
set maxamp2=`echo $maxamp | awk '{ print 2*$1}'`\n\
\n\
set timerange=( `awk 'NR==1 { print $1 } { lastx=$1 } END { print lastx }' ${root}_rad.xy` )\n\
\n\
set psfile=${root}-colour.ps\n\
\n\
gmtdefaults -D >.gmtdefaults4\n\
gmtset PAGE_ORIENTATION portrait MEASURE_UNIT cm WANT_EURO_FONT TRUE LABEL_FONT_SIZE 12 ANOT_FONT_SIZE 10 PAPER_MEDIA a4 D_FORMAT %%lg LABEL_OFFSET 0.05c ANNOT_OFFSET_PRIMARY 0.15c\n\
# 3cm Descriptive text\n\
pstext -M -X0 -Y0 -R0/20/0/30 -Jx1 -K > $psfile <${root}.description\n\
\n\
# 8 cm Error surface\n\
grdcontour -X2 -Y22 ${root}_err.grd -C${root}.cont -R$grdrange -JX9.1/5.8 -B0.5:\"Splitting Delay (s)\":/20:\"Fast direction\":WSen -O -K -A-1f1 -G1000 -Wa1.5p -Wc0.5p >>$psfile\n\
psxy -R -JX -S+0.3 -W3p,150/150/150  -O -K  >>$psfile <<EOF\n\
$besttime $bestfast\n\
EOF\n\
psxy -R -JX -W1p/200/200/200to -O -K  >>$psfile <<EOF\n\
0 $baz\n\
10 $baz\n\
EOF\n\
if ( $?pol ) then \n\
  set polp=`echo $pol | awk '{ if ($1<90) {print $1+90;} else print $1-90 }'`\n\
  psxy -R -JX -M -W1p/200/200/200 -O -K  >>$psfile <<EOF\n\
>\n\
0 $pol\n\
10 $pol\n\
> \n\
0 $polp\n\
10 $polp\n\
EOF\n\
endif\n\
# 4 cm RefRad refTransverse\n\
psxy ${root}_refrad.xy -Y-5.5 -X-1 -R$timerange[1]/$timerange[2]/-$maxamp/$maxamp -JX7.1/3 -Ba5f1/${maxamp}::wseN -W1p -O -K >>$psfile\n\
psxy ${root}_reftra.xy -R$timerange[1]/$timerange[2]/-$maxamp/$maxamp -JX -W1p,blue,- -O -K >>$psfile\n\
pstext <<EOF -JX -R0/1/0/1 -O -K >>$psfile\n\
0.05 0.95 12 0 0 LT Reference Radial-Transverse\n\
EOF\n\
paste ${root}_refrad.xy ${root}_reftra.xy | awk '{ print $4,$2 }' | psxy -X7.1 -R-${maxamp}/${maxamp}/-${maxamp}/${maxamp} -JX3/3 -W1p -B${maxamp}:\"Transverse\":/${maxamp2}:\"Radial\":wsen  -O -K >>$psfile\n\
\n\
# 4 cm Data radial-transverse\n\
psxy ${root}_rad.xy -X-7.1 -Y-3 -R$timerange[1]/$timerange[2]/-$maxamp/$maxamp -JX7.1/3 -Ba5f1/${maxamp}::wsen -W1p -O -K >>$psfile\n\
psxy ${root}_tra.xy -R$timerange[1]/$timerange[2]/-$maxamp/$maxamp -JX -W1p,blue,- -O -K >>$psfile\n\
pstext <<EOF -JX -R0/1/0/1 -O -K >>$psfile\n\
0.05 0.95 12 0 0 LT Observed Radial-Transverse\n\
EOF\n\
paste ${root}_rad.xy ${root}_tra.xy | awk '{ print $4,$2 }' | psxy -X7.1  -R-${maxamp}/${maxamp}/-${maxamp}/${maxamp} -JX3/3 -W1p -B${maxamp}:\"Transverse\":/${maxamp2}:\"Radial\":wsen  -O -K >>$psfile\n\
\n\
# 4 cm Fast Slow\n\
psxy ${root}_fast.xy -X-7.1 -Y-3 -R$timerange[1]/$timerange[2]/-$maxamp/$maxamp -JX7.1/3 -Ba5f1/${maxamp}::wsen -W1p -O -K >>$psfile\n\
psxy ${root}_slow.xy -R$timerange[1]/$timerange[2]/-$maxamp/$maxamp -JX -W1p,blue,- -O -K >>$psfile\n\
pstext <<EOF -JX -R0/1/0/1 -O -K >>$psfile\n\
0.05 0.95 12 0 0 LT Observed Fast-Slow\n\
EOF\n\
paste ${root}_fast.xy ${root}_slow.xy | awk '{ print $4,$2 }' | psxy -X7.1  -R-${maxamp}/${maxamp}/-${maxamp}/${maxamp} -JX3/3 -W1p -B${maxamp}:\"Slow\":/${maxamp2}:\"Fast\":wsen  -O -K >>$psfile\n\
#psxy ${root}_fs_pmp.xy -X7.1  -R-${maxamp}/${maxamp}/-${maxamp}/${maxamp} -JX3/3 -W1p -B${maxamp}:\"Slow\":/${maxamp2}:\"Fast\":wsen  -O -K >>$psfile\n\
\n\
# 4 cm Fast Slow, corrected\n\
psxy ${root}_fastcor.xy -X-7.1 -Y-3 -R$timerange[1]/$timerange[2]/-$maxamp/$maxamp -JX7.1/3 -Ba5f1/${maxamp}::wsen -W1p -O -K >>$psfile\n\
psxy ${root}_slowcor.xy -R$timerange[1]/$timerange[2]/-$maxamp/$maxamp -JX -W1p,blue,- -O -K >>$psfile\n\
pstext <<EOF -JX -R0/1/0/1 -O -K >>$psfile\n\
0.05 0.95 12 0 0 LT Corrected Fast-Slow\n\
EOF\n\
paste ${root}_fastcor.xy ${root}_slowcor.xy | awk '{ print $4,$2 }' | psxy -X7.1  -R-${maxamp}/${maxamp}/-${maxamp}/${maxamp} -JX3/3 -W1p -B${maxamp}:\"Slow\":/${maxamp2}:\"Fast\":wsen  -O -K >>$psfile\n\
#psxy ${root}_fscor_pmp.xy -X7.1  -R-${maxamp}/${maxamp}/-${maxamp}/${maxamp} -JX3/3 -W1p -B${maxamp}:\"Slow\":/${maxamp2}:\"Fast\":wsen  -O -K >>$psfile\n\
\n\
\n\
# 4 cm Rad Transverse, corrected \n\
psxy ${root}_radcor.xy -X-7.1 -Y-3 -R$timerange[1]/$timerange[2]/-$maxamp/$maxamp -JX7.1/3 -Ba5f1:\"Time (s)\":/${maxamp}::wsen -W1p -O -K >>$psfile\n\
psxy ${root}_tracor.xy -R$timerange[1]/$timerange[2]/-$maxamp/$maxamp -JX -W1p,blue,- -O -K >>$psfile\n\
pstext <<EOF -JX -R0/1/0/1 -O -K >>$psfile\n\
0.05 0.95 12 0 0 LT Corrected Radial-Transverse\n\
EOF\n\
paste ${root}_radcor.xy ${root}_tracor.xy | awk '{ print $4,$2 }' | psxy -X7.1 -R-${maxamp}/${maxamp}/-${maxamp}/${maxamp} -JX3/3 -W1p -B${maxamp}:\"Transverse\":/${maxamp2}:\"Radial\":wsen  -O -K >>$psfile\n\
#psxy ${root}_rtcor_pmp.xy -X7.1 -R-${maxamp}/${maxamp}/-${maxamp}/${maxamp} -JX3/3 -W1p -B${maxamp}:\"Transverse\":/${maxamp2}:\"Radial\":wsen  -O -K >>$psfile\n\
\n\
# make amplitude and delay corrected traces\n\
awk '{ print $1-delay,$2*alpha/(1-alpha) }' delay=$delay alpha=$alpha ${root}_radcor.xy > ${root}_tmp_radcorcor.xy\n\
awk '{ print $1-delay,$2*alpha/(1-alpha) }' delay=$delay alpha=$alpha ${root}_tracor.xy > ${root}_tmp_tracorcor.xy\n\
# Interpolate reference station to grid of those data\n\
sample1d  ${root}_refrad.xy -N${root}_tmp_radcorcor.xy > ${root}_tmp_refrad.xy\n\
sample1d  ${root}_reftra.xy -N${root}_tmp_tracorcor.xy > ${root}_tmp_reftra.xy\n\
# Resample corrected data again on this grid to make sure we really share the same grid\n\
sample1d  ${root}_tmp_radcorcor.xy -N${root}_tmp_refrad.xy > ${root}_tmp_radcc.xy\n\
sample1d  ${root}_tmp_tracorcor.xy -N${root}_tmp_reftra.xy > ${root}_tmp_tracc.xy\n\
# Residual Observed-reference\n\
paste ${root}_tmp_radcc.xy ${root}_tmp_refrad.xy | awk '{ print $1,$2-$4 }' > ${root}_tmp_resrad.xy\n\
paste ${root}_tmp_tracc.xy ${root}_tmp_reftra.xy | awk '{ print $1,$2-$4 }' > ${root}_tmp_restra.xy\n\
\n\
\n\
# 4 cm Residual (Observed - reference)\n\
psxy ${root}_tmp_resrad.xy -X-7.1 -Y-3 -R$timerange[1]/$timerange[2]/-$maxamp/$maxamp -JX7.1/3 -Ba5f1:\"Time (s)\":/${maxamp}::wSen -W1p -O -K >>$psfile\n\
psxy ${root}_tmp_restra.xy -R$timerange[1]/$timerange[2]/-$maxamp/$maxamp -JX -W1p,blue,- -O -K >>$psfile\n\
pstext <<EOF -JX -R0/1/0/1 -O -K >>$psfile\n\
0.05 0.95 12 0 0 LT Residual Radial-Transverse\n\
EOF\n\
paste ${root}_tmp_resrad.xy ${root}_tmp_restra.xy | awk '{ print $4,$2 }' | psxy -X7.1 -R-${maxamp}/${maxamp}/-${maxamp}/${maxamp} -JX3/3 -W1p -B${maxamp}:\"Transverse\":/${maxamp2}:\"Radial\":wsen  -O -K >>$psfile\n\
#psxy ${root}_rtcor_pmp.xy -X7.1 -R-${maxamp}/${maxamp}/-${maxamp}/${maxamp} -JX3/3 -W1p -B${maxamp}:\"Transverse\":/${maxamp2}:\"Radial\":wsen  -O -K >>$psfile\n\
\n\
psxy < /dev/null -Jx1 -R -O >>$psfile\n\
\\rm ${root}_tmp*\n\
");
      fclose(output);

      VRB(printf("Executing GMT script\n"));
      sprintf(cmdstring,"csh %s-colour.gmt",par->root);
      status=system(cmdstring);
      if (status) {
	sprintf(warn_str,"Execution of GMT script failed. Error status: %d",status);
	warn_msg(warn_str);
      }

    }

    gsl_vector_float_free(vec1);
    gsl_vector_float_free(vec2);
  }
}

FILE *open_for_write(char *root,char *extension){
  char tmpstring[256];
  FILE *output;
  strcpy(tmpstring, root);
  strcat(tmpstring, extension);
  
  output=fopen(tmpstring,"wb");
  if (!output) { 
    sprintf(abort_str,"Cannot open %s for output.",tmpstring);
    abort_msg(abort_str);
  }
  return(output);
}


void make_event_name(char *string, sachdr *hdr, int mode) {
  /* make event name from information in SAC header.
     mode=EVN_YYJJJHHMM    yy jjj hh:mm */
 switch(mode) {
 case EVN_YYJJJHHMM:
   sprintf(string,"%02d %03d %02d:%02d",hdr->nzyear%100,hdr->nzjday,hdr->nzhour,hdr->nzmin);
   break;
 default:
   abort_msg("make_event_name: Illegal mode");
 }
}


void single_split_correl(int method, hor_split *hsplit, float maxshift, gsl_vector_float *nspl, gsl_vector_float *espl, long beg, gsl_vector_float *ref_north, gsl_vector_float *ref_east, long ref_beg, long len, float delta, gsl_matrix *m_res_energy, gsl_matrix *m_delay, gsl_matrix *m_alpha) {
  /* Grid search for splitting by comparison with reference trace */
  /* transcribed from matlab function split/grdsplit5.m */
  long maxlag,maxxclag,npts,itime,itshiftf,itshifts;
  int i,j,k;
  double fast,time;
  double mc,energy_ref,energy_obs,alpha,en_res;
  long mi,idelay;

  gsl_vector_float_view vue1,vue2,nrefview,erefview;
  gsl_vector_view dvue1,dvue2;
  gsl_vector_float *nref;
  gsl_vector_float *eref;
  gsl_vector *x1;
  gsl_vector *EOns[2][2],*Cns[2][2],*EOfs[2][2],*Cfs[2][2];
  double nsplfi, nsplbi_1, nsplf1_i, nsplb_i, esplfi, esplbi_1, esplf1_i, esplb_i ;


  npts=nspl->size;
  if (beg==0 && len==npts) 
   abort_msg("Fourier-domain calculation of correlation for reference station method not yet implemented");
  

  /* Isolate relevant part of reference vector from ref_bef Length len*/
  nrefview=gsl_vector_float_subvector(ref_north,ref_beg,len);
  erefview=gsl_vector_float_subvector(ref_east,ref_beg,len);
  nref=&nrefview.vector;
  eref=&erefview.vector;

  /* xc=zeros(length(vfstdir),length(vtsplt)) ; */
  /* delay=zeros(size(xc)); */
  /* en_res=zeros(size(xc)); */

  /* maxlag=round(MAXLAG/SAMPLE_INTERVAL); */
  maxlag=(long)ROUND(maxshift/delta);

  /* maxxclag=maxlag+ceil(0.5*max(vtsplt)/SAMPLE_INTERVAL) ; */
  maxxclag=maxlag+(long)ROUND(hsplit->timemax/delta) ;

  VRB(printf("maxlag %ld maxxclag %ld\n",maxlag,maxxclag));
  /* len=length(nref) ; */
  /* f=bspl+len-1 ; */
  /* overlag=max(maxxclag-bspl+1,bspl+len-1+maxxclag-length(nspl)) ; */
  /* if overlag>0 */
  /*   maxlag=maxlag-overlag ; */
  /*   maxxclag=maxxclag-overlag ; */
  /*   if(maxxclag<0) */
  /*      error('nspl,espl too short to accommodate maximum time shift even at zero delay') */
  /*   else */
  /*      warning(sprintf('nspl,espl too short to accommodate maximum time shift and lag without introducing edge effects\nReducing maxlag to %4.2f sec.',maxlag*SAMPLE_INTERVAL)) */
  /*   end */
  /* end */

  x1=gsl_vector_alloc(2*maxlag+1);

  /* energy_ref=nref'*nref + eref'*eref ;  */
  energy_ref=gsl_df_normsqr(nref)+gsl_df_normsqr(eref);
  VRB(printf("Energy Ref: %g\n",energy_ref));
  /*   % energy_obs has to be calculated for all possible lags (implying a different window) and separately for n,e as these are potentially time-shifted */
  /* % Matrix of energy in observed signal for different window shift, have to calculate n*n,e*e, and n*e to be able to */
  /* % transform into fast,slow direction */
  /* EOns=zeros(2*maxxclag+1,2,2); */
  /*   EOns[0][0]=gsl_vector_alloc(2*maxxclag+1); */
  /*   EOns[0][1]=gsl_vector_alloc(2*maxxclag+1); */
  /*   EOns[1][0]=gsl_vector_alloc(2*maxxclag+1); */
  /*   EOns[1][1]=gsl_vector_alloc(2*maxxclag+1); */
  /*   EOfs[0][0]=gsl_vector_alloc(2*maxxclag+1); */
  /*   EOfs[0][1]=gsl_vector_alloc(2*maxxclag+1); */
  /*   EOfs[1][0]=gsl_vector_alloc(2*maxxclag+1); */
  /*   EOfs[1][1]=gsl_vector_alloc(2*maxxclag+1); */
  /*   Cns[0][0] =gsl_vector_alloc(2*maxxclag+1); */
  /*   Cns[0][1] =gsl_vector_alloc(2*maxxclag+1); */
  /*   Cns[1][0] =gsl_vector_alloc(2*maxxclag+1); */
  /*   Cns[1][1] =gsl_vector_alloc(2*maxxclag+1); */
  /*   Cfs[0][0] =gsl_vector_alloc(2*maxxclag+1); */
  /*   Cfs[0][1] =gsl_vector_alloc(2*maxxclag+1); */
  /*   Cfs[1][0] =gsl_vector_alloc(2*maxxclag+1); */
  /*   Cfs[1][1] =gsl_vector_alloc(2*maxxclag+1); */
  CMAT_ALLOCATE(EOns,2*maxxclag+1);
  CMAT_ALLOCATE(EOfs,2*maxxclag+1);
  CMAT_ALLOCATE(Cns,2*maxxclag+1);
  CMAT_ALLOCATE(Cfs,2*maxxclag+1);

  if (beg==0 && len==npts) {
    abort_msg("FFT type cross-correlation for reference station correlation method not yet implemented");
  } else {
    /* EOns(maxxclag+1,1,1)=nspl(bspl:f)'*nspl(bspl:f) ; */
    vue1=gsl_vector_float_subvector(nspl,beg,len);
    vue2=gsl_vector_float_subvector(espl,beg,len);
    gsl_vector_set(EOns[0][0],maxxclag,gsl_df_normsqr(&vue1.vector));
    /* EOns(maxxclag+1,2,2)=espl(bspl:f)'*espl(bspl:f) ; */
    gsl_vector_set(EOns[1][1],maxxclag,gsl_df_normsqr(&vue2.vector));
    /* EOns(maxxclag+1,2,1)=nspl(bspl:f)'*espl(bspl:f) ; */
    gsl_vector_set(EOns[1][0],maxxclag,gsl_df_dotprod(&vue1.vector,&vue2.vector));
    gsl_vector_set(EOns[0][1],maxxclag,gsl_df_dotprod(&vue1.vector,&vue2.vector));
    /* determine correlations for other window positions by sequential updates */
    VRB(printf("EOns(zero lag): [ %g %g ; %g %g ]\n",
	       gsl_vector_get(EOns[0][0],maxxclag),
	       gsl_vector_get(EOns[0][1],maxxclag),
	       gsl_vector_get(EOns[1][0],maxxclag),
	       gsl_vector_get(EOns[1][1],maxxclag)));
    /* for i=1:maxxclag */
    for(i=1;i<=maxxclag;i++){
      nsplfi  =(double)gsl_vector_float_get(nspl,beg+len-1+i); /* nspl(f+i) */
      nsplbi_1=(double)gsl_vector_float_get(nspl,beg+i-1);     /* nspl(bspl+i-1) */
      nsplf1_i=(double)gsl_vector_float_get(nspl,beg+len-i);   /* nspl(f+1-i) */
      nsplb_i =(double)gsl_vector_float_get(nspl,beg-i);       /* nspl(bspl-i) */
      esplfi  =(double)gsl_vector_float_get(espl,beg+len-1+i); /* espl(f+i) */
      esplbi_1=(double)gsl_vector_float_get(espl,beg+i-1);     /* espl(bspl+i-1) */
      esplf1_i=(double)gsl_vector_float_get(espl,beg+len-i);   /* espl(f+1-i) */
      esplb_i =(double)gsl_vector_float_get(espl,beg-i);       /* espl(bspl-i) */
  
      /*   EOns(maxxclag+1+i,1,1)=EOns(maxxclag+i,1,1)  +nspl(f+i)^2  -nspl(bspl+i-1)^2 ; */
      gsl_vector_set(EOns[0][0],maxxclag+i,
		     gsl_vector_get(EOns[0][0],maxxclag+i-1)+SQR(nsplfi)-SQR(nsplbi_1));
      /*   EOns(maxxclag+1-i,1,1)=EOns(maxxclag+2-i,1,1)-nspl(f+1-i)^2+nspl(bspl-i)^2   ; */
      gsl_vector_set(EOns[0][0],maxxclag-i,
		     gsl_vector_get(EOns[0][0],maxxclag-i+1)-SQR(nsplf1_i)+SQR(nsplb_i));
      /*   EOns(maxxclag+1+i,2,2)=EOns(maxxclag+i,2,2)  +espl(f+i)^2  -espl(bspl+i-1)^2 ; */
      gsl_vector_set(EOns[1][1],maxxclag+i,
		     gsl_vector_get(EOns[1][1],maxxclag+i-1)+SQR(esplfi)-SQR(esplbi_1));
      /*   EOns(maxxclag+1-i,2,2)=EOns(maxxclag+2-i,2,2)-espl(f+1-i)^2+espl(bspl-i)^2   ; */
      gsl_vector_set(EOns[1][1],maxxclag-i,
		     gsl_vector_get(EOns[1][1],maxxclag-i+1)-SQR(esplf1_i)+SQR(esplb_i));
      /*   EOns(maxxclag+1+i,2,1)=EOns(maxxclag+i,2,1)  +nspl(f+i)*espl(f+i)  -nspl(bspl+i-1)*espl(bspl+i-1) ; */
      gsl_vector_set(EOns[1][0],maxxclag+i,
		     gsl_vector_get(EOns[1][0],maxxclag+i-1)+nsplfi*esplfi-nsplbi_1*esplbi_1);
      /*   EOns(maxxclag+1-i,2,1)=EOns(maxxclag+2-i,2,1)-nspl(f+1-i)*espl(f+1-i)+nspl(bspl-i)*espl(bspl-i) ; */
      gsl_vector_set(EOns[1][0],maxxclag-i,
		     gsl_vector_get(EOns[1][0],maxxclag-i+1)-nsplf1_i*esplf1_i+nsplb_i*esplb_i);
      /* end */
      /* EOns(:,1,2)=EOns(:,2,1) ; */
      gsl_vector_set(EOns[0][1],maxxclag+i,gsl_vector_get(EOns[1][0],maxxclag+i));
      gsl_vector_set(EOns[0][1],maxxclag-i,gsl_vector_get(EOns[1][0],maxxclag-i));
    }
    /* plot out EOns */

    /* % calculate elements of cross-correlation matrix */
    /* Cns=zeros(2*maxxclag+1,2,2); */
    Cns[0][0]=gsl_vector_alloc(2*maxxclag+1);
    Cns[0][1]=gsl_vector_alloc(2*maxxclag+1);
    Cns[1][0]=gsl_vector_alloc(2*maxxclag+1);
    Cns[1][1]=gsl_vector_alloc(2*maxxclag+1);
    /* Cns(:,1,1)=corrshift(nref,nspl,bspl,maxxclag);  */
    gsl_df_corrshift(Cns[0][0],nref,0,nspl,beg,len,maxxclag,GSL_CORRSHIFT_TWO);
    /* Cns(:,1,2)=corrshift(eref,nspl,bspl,maxxclag); */
    gsl_df_corrshift(Cns[0][1],eref,0,nspl,beg,len,maxxclag,GSL_CORRSHIFT_TWO);
    /* Cns(:,2,1)=corrshift(nref,espl,bspl,maxxclag);  */
    gsl_df_corrshift(Cns[1][0],nref,0,espl,beg,len,maxxclag,GSL_CORRSHIFT_TWO);
    /* Cns(:,2,2)=corrshift(eref,espl,bspl,maxxclag); */
    gsl_df_corrshift(Cns[1][1],eref,0,espl,beg,len,maxxclag,GSL_CORRSHIFT_TWO);

/*     {  */
/*       FILE *tmpfid; */
/*       tmpfid=fopen("tmp.Cns00","w"); */
/*       gsl_vector_fprintf(tmpfid, Cns[0][0],"%g"); */
/*       fclose(tmpfid); */
/*       tmpfid=fopen("tmp.Cns01","w"); */
/*       gsl_vector_fprintf(tmpfid, Cns[0][1],"%g"); */
/*       fclose(tmpfid); */
/*       tmpfid=fopen("tmp.Cns10","w"); */
/*       gsl_vector_fprintf(tmpfid, Cns[1][0],"%g"); */
/*       fclose(tmpfid); */
/*       tmpfid=fopen("tmp.Cns11","w"); */
/*       gsl_vector_fprintf(tmpfid, Cns[1][1],"%g"); */
/*       fclose(tmpfid); */
/*     } */

  }  
  /* disp('Start grid search'); */
  /* for i=1:length(vfstdir) */
  /*   fstdir=vfstdir(i) ; */
  for (fast=hsplit->fastmin,j=0; fast<=hsplit->fastmax+TOLERANCE; fast+=hsplit->faststep, j++) {
    VRB(printf("j=%d Fast=%f \n",j,fast));
    /*   Cfs=rotcr(Cns,fstdir); */
/*     gsl_vector_memcpy(Cfs[0][0],Cns[0][0]); */
/*     gsl_vector_memcpy(Cfs[0][1],Cns[0][1]); */
/*     gsl_vector_memcpy(Cfs[1][0],Cns[1][0]); */
/*     gsl_vector_memcpy(Cfs[1][1],Cns[1][1]); */
    CMAT_MEMCPY(Cfs,Cns);
    gsl_rotcr(Cfs,fast);
    /* % the following statement is a bit of an overkill as only few elements of EOns are actually needed. However, */
    /* % it makes the code look somewhat neater and will actually be faster   */
    /* % if length(vtsplt) is of similar magnitude to maxxclag+length(nref) */
    /*   EOfs=rotcr(EOns,fstdir); */
    CMAT_MEMCPY(EOfs,EOns);
    gsl_rotcr(EOfs,fast);

/*   for j=1:length(vtsplt) */
/*     tsplt=vtsplt(j) ; */
    for (time=hsplit->timemin,k=0; time<=hsplit->timemax+TOLERANCE; time+=hsplit->timestep, k++) {
      /*       tshift=round(tsplt/(2*SAMPLE_INTERVAL)) ; */
      itime=(long)ROUND(time/delta);
/*       VRB(printf("k=%d Time=%f itime=%d\n",k,time,itime)); */
      /*     if tsplt==0.0 & i>=2 */
      /* %     for zero splitting, fast direction does not matter and we can */
      /* %     just copy the result from the first fast direction */
      if(itime==0 && j>=1) {
	/* fast direction is irrelevant for zero splitting time, ie we can just copy result */
	/*       en_res(i,j)=en_res(1,j) ; */
	/* m_res_energy[j,k]=m_res_energy[1,k] */
	gsl_matrix_set(m_res_energy,j,k,gsl_matrix_get(m_res_energy,0,k));
/*       xc(i,j)=xc(1,j) ;  */

/*       delay(i,j)=delay(1,j) ; */
	gsl_matrix_set(m_delay,j,k,gsl_matrix_get(m_delay,0,k));

	gsl_matrix_set(m_alpha,j,k,gsl_matrix_get(m_alpha,0,k));
	continue;
      /*     else */
      } else {
	/*       tshift=round(tsplt/(2*SAMPLE_INTERVAL)) ; */
       	/*       x1=tshiftadd(Cfs(:,1,1),Cfs(:,2,2),tshift,maxlag); */
	itshiftf=itime/2;      
	itshifts=itime-itshiftf;
	dvue1=gsl_vector_subvector(Cfs[0][0],maxxclag-maxlag-itshiftf,2*maxlag+1);
	dvue2=gsl_vector_subvector(Cfs[1][1],maxxclag-maxlag+itshifts,2*maxlag+1);
	gsl_vector_memcpy(x1,&dvue1.vector);
	gsl_vector_add(x1,&dvue2.vector);
	/*       [mc,mi]=max(x1) ; */
	mi=(long)gsl_vector_max_index(x1);
	mc=gsl_vector_get(x1,mi);
	/* %%%      xc(i,j)=mc ; */
	/*       delay(i,j)=(mi-(length(x1)+1)/2); */
	idelay=(mi-maxlag);
	gsl_matrix_set(m_delay,j,k,idelay*delta);
	/* % energy_obs: energy in spl seismogram within reference window after applying delay and splitting correction */
	/*       me=maxxclag+1+delay(i,j); */
	/*       energy_obs=EOfs(me-tshift,1,1)+EOfs(me+tshift,2,2) ; */
	energy_obs=gsl_vector_get(EOfs[0][0],maxxclag-itshiftf+idelay)
                  +gsl_vector_get(EOfs[1][1],maxxclag+itshiftf+idelay);
/* 	%    calculate residual energy between corrected seismogram and reference seismogram allowing for amplitude differences (final amplitude is normalised with respect to amplitude of reference signal */
/* %     NB mc(fast,slow)=mc(north,east) */
/* %        can calculate residual energy directly, without */
/* %        calculating alpha (amplitude correction factor) and residual */
/* %        vectors first, after some algebra (in notes). Prefactor 4 such that this corresponds to physical */
/* %        energy if there is no amplitude difference */
/*       en_res(i,j)=4*(energy_ref*energy_obs^2+2*energy_ref*energy_obs*mc-(energy_ref+energy_obs)*mc^2 ... */
/* 	         -2*mc^3+energy_ref^2*energy_obs)/(energy_ref+2*mc+energy_obs)^2; */
	en_res=4*(energy_ref*SQR(energy_obs)+2*energy_ref*energy_obs*mc-(energy_ref+energy_obs)*SQR(mc)-2*CUB(mc)+SQR(energy_ref)*energy_obs)/SQR(energy_ref+2*mc+energy_obs);

/* % ALTERNATIVE : no alpha correction */
/* %      en_res(i,j)=(energy_ref-2*mc+energy_obs); */
/*       if(lshow) */
/* %        Axis scale factors are adjusted such that reference, corrected data, and residual time series are  */
/* %          visually directly comparable, i.e. amplitude correction has been performed. They are in a way */
/* %          normalised to amplitude of reference signal (but not actually, just on figure) */
/*          scale=1.1*max(abs([nref ; eref])); */
/*          [ns,es]=splitopdel(nspl,espl,bspl,len,tsplt,fstdir,delay(i,j)*SAMPLE_INTERVAL,-1) ; */
/*          [pol,lin]=pmotion(ns,es) ; */

/*          alpha=(energy_ref+mc)/(energy_obs+2*mc+energy_ref) ;  */
      alpha=(energy_ref+mc)/(energy_obs+2*mc+energy_ref);
      /* in contrast to matlab version, normalise by total energy available */
      gsl_matrix_set(m_res_energy,j,k,en_res/(SQR(1-alpha)*energy_ref+SQR(alpha)*energy_obs));
/*       gsl_matrix_set(m_res_energy,j,k,en_res); */
      gsl_matrix_set(m_alpha,j,k,alpha) ;

/* % ALTERNATIVE : no alpha correction */
/* %	alpha=0.5 ;  */
/* %       (alpha < 0.5   amplitudes of reference trace smaller than station trace */
/* %       (alpha > 0.5   amplitudes of reference trace larger than station race */
/*          nres=2*((1-alpha)*nref-alpha*ns); */
/*          eres=2*((1-alpha)*eref-alpha*es); */
/*          fprintf('DEBUG: Direct obs. energy estimate :%.4g vs. pre-calc %.4g dif = %.4g\n', ... */
/* 	      nres'*nres+eres'*eres,en_res(i,j),nres'*nres+eres'*eres-en_res(i,j)) ; */
/*          fprintf('DEBUG: Direct obs.energy :%.4g vs. precalc %.4g dif = %.4g\n', ... */
/*               ns'*ns+es'*es,energy_obs,ns'*ns+es'*es-energy_obs) ; */
/*          fprintf('DEBUG: Direct cross.correl :%.4g vs. precalc %.4g dif = %.4g\n', ... */
/* 	      nref'*ns+eref'*es,mc,nref'*ns+eref'*es-mc) ; */
/* 	 if(en_res(i,j)<0) */
/* 	   error('Negative energy. This indicates a bug in the program') */
/* 	 end */
	VRB(if (k==1)printf(" m_res=%g     delay=%f  alpha=%f\n",gsl_matrix_get(m_res_energy,j,k),gsl_matrix_get(m_delay,j,k),gsl_matrix_get(m_alpha,j,k)));
      }
    } /* end loop over splitting times */
  } /* end loop over fast directions */
  /*   gsl_vector_free(EOns[0][0]); */
  /*   gsl_vector_free(EOns[0][1]); */
  /*   gsl_vector_free(EOns[1][0]); */
  /*   gsl_vector_free(EOns[1][1]); */
  /*   gsl_vector_free(EOfs[0][0]); */
  /*   gsl_vector_free(EOfs[0][1]); */
  /*   gsl_vector_free(EOfs[1][0]); */
  /*   gsl_vector_free(EOfs[1][1]); */
  /*   gsl_vector_free(Cns[0][0] ); */
  /*   gsl_vector_free(Cns[0][1] ); */
  /*   gsl_vector_free(Cns[1][0] ); */
  /*   gsl_vector_free(Cns[1][1] ); */
  /*   gsl_vector_free(Cfs[0][0] ); */
  /*   gsl_vector_free(Cfs[0][1] ); */
  /*   gsl_vector_free(Cfs[1][0] ); */
  /*   gsl_vector_free(Cfs[1][1] ); */
  CMAT_FREE(EOns);
  CMAT_FREE(EOfs);
  CMAT_FREE(Cns);
  CMAT_FREE(Cfs);
}

void single_split_sks(int method, hor_split *hsplit, gsl_vector_float *north, gsl_vector_float *east, long beg, long len, float delta, float baz, gsl_matrix *m_res_energy, gsl_matrix *m_pol) {
 /* Grid search for SKS splitting */
  float *cnn,*cee,*cne,*data_e,*data_n;
  long npts;
  long maxlag,n2,itime;
  int i,j,k;
  double fast,c,s,time;
  double pol,pol_n,pol_e,pol_f,pol_s,tot_energy;
  /* for matrix calculations */
  gsl_matrix *m2_cor0=gsl_matrix_alloc(2,2);
  gsl_matrix *m2_corn=gsl_matrix_alloc(2,2);
  gsl_matrix *m2_rot=gsl_matrix_alloc(2,2);
  gsl_matrix *m2_rotbaz=gsl_matrix_alloc(2,2);
  gsl_matrix *m2_aux=gsl_matrix_alloc(2,2);
  gsl_matrix *m2_fscor0=gsl_matrix_alloc(2,2);
  gsl_matrix *m2_fscorn=gsl_matrix_alloc(2,2);
  /* for eigenvalue analysis */
  gsl_eigen_symmv_workspace *ws_eigen=gsl_eigen_symmv_alloc(2);
  gsl_matrix *m2_evec=gsl_matrix_alloc(2,2);
  gsl_vector *eval=gsl_vector_alloc(2);

  VRB(printf("in single_split_sks %f %f %f\n",hsplit->timemax,delta,ROUND(hsplit->timemax/delta)));
  npts=north->size;
  maxlag=(long)ROUND(hsplit->timemax/delta);

  /* Calculate cross-correlation matrices */
  if (beg==0 && len==npts) {
/*     FILE *tmpfid; */
/*     gsl_vector_float_view tmpvecview; */
    /* tapered sequences - can use fft to calculate cross-correlations */
    n2=nxtpwr2(npts+maxlag);
    VRB(printf("calculate cross-correlation with FFT maxlag=%ld npts=%ld n2=%ld\n",maxlag,npts,n2));
    cnn=(float *)malloc(n2*sizeof(float));
    cee=(float *)malloc(n2*sizeof(float));
    cne=(float *)malloc(n2*sizeof(float));

    data_n=(float *)calloc(sizeof(float),n2); 
    data_e=(float *)calloc(sizeof(float),n2);
    if (! ( cnn || cee || cne || data_n || data_e )) {
      abort_msg("single-split-sks: Out of memory");
    }
    memcpy(data_n,north->data,npts*sizeof(float));
    memcpy(data_e,east->data,npts*sizeof(float));
    /* FFT */
/*     VRB(printf("FFT\n")); */
    gsl_fft_real_float_radix2_transform(data_n,1,n2);
    gsl_fft_real_float_radix2_transform(data_e,1,n2);
    /* Multiply with complex conjugate in Freq domain */
/*     VRB(printf("Correlation in freq domain\n")); */
    for (i=1;i<=n2/2;i++) {
      cnn[i]=ABS2(data_n[i],data_n[n2-i]); /*REAL*/   /* cnn=data_n* data_n */
      cnn[n2-i]=0;              /*IMAG*/
      cee[i]=ABS2(data_e[i],data_e[n2-i]); /*REAL*/   /* cee=data_e* data_e */
      cee[n2-i]=0;             /*IMAG*/
      cne[i]=data_n[i]*data_e[i]+data_n[n2-i]*data_e[n2-i]; /*REAL*/  /* cne=data_n* data_e */
      cne[n2-i]=-data_n[n2-i]*data_e[i]+data_n[i]*data_e[n2-i];              /*IMAG*/
    }
    /* IFFT */
/*     VRB(printf("IFFT\n");) */
    gsl_fft_halfcomplex_float_radix2_inverse(cnn,1,n2);
    gsl_fft_halfcomplex_float_radix2_inverse(cee,1,n2);
    gsl_fft_halfcomplex_float_radix2_inverse(cne,1,n2);
    /* NB time domain xcorrel stored in wrap around order - only good up to a lag of maxlag */
    /* DEBUG: output cross correlations */
/*     tmpfid=fopen("tmp.cnn","w"); */
/*     tmpvecview=gsl_vector_float_view_array(cnn,n2); */
/*     gsl_vector_float_fprintf(tmpfid, &tmpvecview.vector,"%g"); */
/*     fclose(tmpfid); */

/*     tmpfid=fopen("tmp.cee","w"); */
/*     tmpvecview=gsl_vector_float_view_array(cee,n2); */
/*     gsl_vector_float_fprintf(tmpfid, &tmpvecview.vector,"%g"); */
/*     fclose(tmpfid); */

/*     tmpfid=fopen("tmp.cne","w"); */
/*     tmpvecview=gsl_vector_float_view_array(cne,n2); */
/*     gsl_vector_float_fprintf(tmpfid, &tmpvecview.vector,"%g"); */
/*     fclose(tmpfid); */
    /* END DEBUG */

    tot_energy=cnn[0]+cee[0];   /* normalise by sum of auto-correlation at zero lag, representative of the total energy in the two components */
  }
  else {
      /* calculate by time-shifting window */
    abort_msg("Time-shifting algorithm for SKS not yet implemented");
  }
  gsl_matrix_set(m2_cor0,0,0,(double)cnn[0]);
  gsl_matrix_set(m2_cor0,0,1,(double)cne[0]);
  gsl_matrix_set(m2_cor0,1,0,(double)cne[0]);
  gsl_matrix_set(m2_cor0,1,1,(double)cee[0]);
/*   VRB(printf("m2_cor0:\n"); gsl_matrix_fprintf(stdout,m2_cor0,"%g")); */

  for (fast=hsplit->fastmin,j=0; fast<=hsplit->fastmax+TOLERANCE; fast+=hsplit->faststep, j++) {
    VRB(printf("j=%d Fast=%f \n",j,fast));
    /* rot: rotate from ne to fs */
    c=cos(fast*PI/180); s=sin(fast*PI/180);
    gsl_matrix_set(m2_rot,0,0,(double)c);
    gsl_matrix_set(m2_rot,0,1,(double)s);
    gsl_matrix_set(m2_rot,1,0,(double)-s);
    gsl_matrix_set(m2_rot,1,1,(double)c);
    /* rotbaz: rotate from fs to rt */
    c=cos((baz-fast)*PI/180); s=sin((baz-fast)*PI/180);	
    gsl_matrix_set(m2_rotbaz,0,0,(double)c);
    gsl_matrix_set(m2_rotbaz,0,1,(double)s);
    gsl_matrix_set(m2_rotbaz,1,0,(double)-s);
    gsl_matrix_set(m2_rotbaz,1,1,(double)c);
    /* fscor0 = rot*cor0*rot'  (zero lag cross-correlation matrix in fast-slow coordinate system */
/*     VRB(printf("m2_rot:\n"); gsl_matrix_fprintf(stdout,m2_rot,"%g")); */
/*     VRB(printf("m2_rotbaz:\n"); gsl_matrix_fprintf(stdout,m2_rotbaz,"%g")); */
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,m2_rot,m2_cor0,0.0,m2_aux);
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,m2_aux,m2_rot,0.0,m2_fscor0);
/*     VRB(printf("m2_fscor0:\n"); gsl_matrix_fprintf(stdout,m2_fscor0,"%g")); */
    for (time=hsplit->timemin,k=0; time<=hsplit->timemax+TOLERANCE; time+=hsplit->timestep, k++) {
      itime=(long)ROUND(time/delta);
/*       VRB(printf("k=%d time=%f itime=%d\n",k,time,itime)); */
      if(itime==0 && j>=1) {
	/* fast direction is irrelevant for zero splitting time, ie we can just copy result */
	/* m_res_energy[j,k]=m_res_energy[1,k] */
	gsl_matrix_set(m_res_energy,j,k,gsl_matrix_get(m_res_energy,0,k));
	gsl_matrix_set(m_pol,j,k,gsl_matrix_get(m_pol,0,k));
	continue;
      } else if (itime==0) {
	/* for zero splitting time can copy correlation matrix for zero lag */
	gsl_matrix_memcpy(m2_fscorn,m2_fscor0); 
      }	else {
	/* corn: cross-correlation matrix at n=itime */
	gsl_matrix_set(m2_corn,0,0,(double)cnn[itime]);
	gsl_matrix_set(m2_corn,0,1,(double)cne[itime]);
	gsl_matrix_set(m2_corn,1,0,(double)cne[n2-itime]); /* cne[itime]=cen[-itime]=cen[n2-itime] */
	gsl_matrix_set(m2_corn,1,1,(double)cee[itime]);
/* 	VRB(printf("m2_corn:\n"); gsl_matrix_fprintf(stdout,m2_corn,"%g")); */
	/* fscorn = rot*corn*rot'  (correlation matrix in fast-slow coordinate system between advanced and retarded traces) */
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,m2_rot,m2_corn,0.0,m2_aux);
	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,m2_aux,m2_rot,0.0,m2_fscorn);
/* 	VRB(printf("m2_fscorn:\n"); gsl_matrix_fprintf(stdout,m2_fscorn,"%g")); */
      }
      /* aux=[fscor0[11],fscorn[12];fscorn[12],fscor0[22];  note that fscorn[12] is used twice rather than fscorn[21] as we need f's but not fs' (with f' being the time-advanced fast component and s the slow component */
      gsl_matrix_set(m2_aux,0,0,gsl_matrix_get(m2_fscor0,0,0));
      gsl_matrix_set(m2_aux,0,1,gsl_matrix_get(m2_fscorn,0,1));
      gsl_matrix_set(m2_aux,1,0,gsl_matrix_get(m2_fscorn,0,1));
      gsl_matrix_set(m2_aux,1,1,gsl_matrix_get(m2_fscor0,1,1));
/*       VRB(printf("m2_aux(composite):\n"); gsl_matrix_fprintf(stdout,m2_aux,"%g")); */

      switch(method){
      case MINEVALUE:
	gsl_eigen_symmv(m2_aux, eval, m2_evec, ws_eigen);              /* Calculate eigenvalues and vectors */
	gsl_eigen_symmv_sort(eval, m2_evec, GSL_EIGEN_SORT_VAL_ASC); /* Sort them */

	gsl_matrix_set(m_res_energy,j,k,gsl_vector_get(eval,0)/tot_energy);   /* Remember smallest eigenvalue */
/* 	VRB(printf("Min. E Value energy: %f\n",gsl_vector_get(eval,0)/tot_energy)); */
	/* Rotate eigenvector associated with larger eigenvalue back into N/S system */
	c=cos((-fast)*PI/180); s=sin((-fast)*PI/180);	
	pol_f=gsl_matrix_get(m2_evec,0,1); pol_s=gsl_matrix_get(m2_evec,1,1);
	pol_n= c*pol_f + s*pol_s; 
	pol_e=-s*pol_f + c*pol_s;
	pol=atan2(pol_e,pol_n)*180/PI; if (pol<0) pol+=360;
	gsl_matrix_set(m_pol,j,k,pol);
	break;
      case MINTRANSVERSE:
	/* aux = rotbaz*aux*rotbaz' (using fscorn for intermediate result)  */
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,m2_rotbaz,m2_aux,0.0,m2_fscorn);
	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,m2_fscorn,m2_rotbaz,0.0,m2_aux);
	gsl_matrix_set(m_res_energy,j,k,gsl_matrix_get(m2_aux,1,1)/tot_energy);
	gsl_matrix_set(m_pol,j,k,fmod(baz,180.));
/* 	VRB(printf("Norm. Transverse energy: %f %f\n",gsl_matrix_get(m2_aux,1,1)/tot_energy,gsl_matrix_get(m_res_energy,j,k))); */
	break;
      }
    } /* continue for(time...) */
  } /* continue for(fast...) */
  /* free allocated memory */
  gsl_matrix_free(m2_cor0);
  gsl_matrix_free(m2_corn);
  gsl_matrix_free(m2_rot);
  gsl_matrix_free(m2_rotbaz);
  gsl_matrix_free(m2_aux);
  gsl_matrix_free(m2_fscor0);
  gsl_matrix_free(m2_fscorn);
  gsl_eigen_symmv_free(ws_eigen);
  gsl_matrix_free(m2_evec);
  gsl_vector_free(eval);
  free(cnn); free(cee); free(cne);
  free(data_n); free(data_e);
}



gsl_vector_float *find_window(sachdr *hdr, gsl_vector_float *data, 
		 analysis_window *win, float maxlag, long *beg, long *len,char *phase) {
  /* determines begin and length of analysis window in terms of samples
     carries out tapering on the actual data if requested
     warns if no tapering is employed and there are less than maxlag seconds
     available (necessary for shift-correlation)

     Input:
     hdr, data  : SAC header, data
     win          : Analysis window structure (phase_start,phase_end, offset_start, offset_end, taper)
     maxlag       : maximum amount of lag (in s) to be expected from later analysis

     Output:
     beg,len      : Begin and length of analysis window in samples (0-based)
                    Note that samples outside the analysis window can be important for
		    non-zero lag if tapering is not employed 
    phase:          Decription of phase used

     Returns:
                    pointer to gsl_vector_float
		 

     Depends On: find_phase
  */
  float del=hdr->delta,b=hdr->b;
  float wb,we;
  char phase_dum[9];
  gsl_vector_float_view vec_vue;
  gsl_vector_float *vec=malloc(sizeof(gsl_vector_float));

  wb=find_phase(hdr,win->phase_start,phase)+win->offset_start;
  if (isnan(wb)) {
    fprintf(stderr,"Could not find phase name %s.",win->phase_start); abort_msg(""); }
  we=find_phase(hdr,win->phase_end,phase_dum)+win->offset_end;
  if (isnan(we)) {
    fprintf(stderr,"Could not find phase name %s.",win->phase_end); abort_msg(""); }

  /*  fprintf(stderr,"DEBUG find_window wb b del refwbeg %d ldum1 %d refwln %d ldum2 %d\n",refwbeg,ldum1,refwlen,ldum2); */
  *beg=(long)ROUND((wb-b)/del);
  *len=(long)ROUND((we-wb)/del)+1;
  if ( *beg<0 || *beg+*len > hdr->npts || *len<3 ) {
    sprintf(warn_str,"Illegal win definition  - beg(s)=%f  end(s)=%f  beg(smpl)=%ld length(smpl)=%ld b(sac)=%f delta=%f npts=%d",
	   wb,we,*beg,*len,b,del,hdr->npts);
    abort_msg(warn_str);
  }

  if (win->taper>=0) {
    VRB(fprintf(stderr,"Tapering with taper length %f.\n", win->taper));
    /* Create sub-vector */
    vec_vue=gsl_vector_float_subvector(data,*beg,*len);
    rmean_and_taper(gsl_vector_float_ptr(&vec_vue.vector,0),*len,del,0,win->taper);
    *beg=0;
    *vec=vec_vue.vector;   /* copy actual vector structure to prevent it going out of scope on exiting */
    return(vec);
  } else {
    /* Check that there is enough data around analysis window to accommodate maximum shifts */
    if (wb-b < maxlag || hdr->e-we < maxlag ) {
      sprintf(warn_str,"Length of time series too short to accommodate maximum lag %f (pre-window: %f, post-window: %f)",
	      maxlag, wb-b, hdr->e-we);
      warn_msg(warn_str);
    }
    return(data);
  }
}


float find_phase(sachdr *hdr,char *phase_name,char *phase_out) {
  /* 
     Input
     hdr:           pointer to SAC header 
     phase_name:    either SAC header variable or phase name (to be matched with pick descriptions
     returns time in s of SAC phase phase_name 
     returns NaN if no match can be found
     Auxiliary output:
     phase_out:         Name of phase finally chosen
   */
  phase_out[8]='\0';
  if (!strcasecmp(phase_name,"b")) {
    strcpy(phase_out," ");
    return(hdr->b); }
  else if (!strcasecmp(phase_name,"e")) {
    strcpy(phase_out," ");
    return(hdr->e); }
  else if (!strcasecmp(phase_name,"o")) {
    strcpy(phase_out,"O");
    return(hdr->o); }
  else if (!strcasecmp(phase_name,"a")) {
/*     if (strncasecmp("-12345",hdr->ka,6))  */
/*       strncpy(phase_out,hdr->ka,8); */
/*     else */
    strcpy(phase_out,"A");
    return(hdr->a); }
  else if (!strcasecmp(phase_name,"f")) {
    strcpy(phase_out,"F");
    return(hdr->f); }
  else if (!strcasecmp(phase_name,"t0")) {
    if (strncasecmp("-12345",hdr->kt0,6)) 
      strncpy(phase_out,hdr->kt0,8);
    else
      strcpy(phase_out,"T0");
    return(hdr->t0); }
  else if (!strcasecmp(phase_name,"t1")) {
    if (strncasecmp("-12345",hdr->kt1,6)) 
      strncpy(phase_out,hdr->kt1,8);
    else
      strcpy(phase_out,"T1");
    return(hdr->t1); }
  else if (!strcasecmp(phase_name,"t2")) {
    if (strncasecmp("-12345",hdr->kt2,6)) 
      strncpy(phase_out,hdr->kt2,8);
    else
      strcpy(phase_out,"T2");
    return(hdr->t2); }
  else if (!strcasecmp(phase_name,"t3")) {
    if (strncasecmp("-12345",hdr->kt3,6)) 
      strncpy(phase_out,hdr->kt3,8);
    else
      strcpy(phase_out,"T3");
    return(hdr->t3); }
  else if (!strcasecmp(phase_name,"t4")) {
    if (strncasecmp("-12345",hdr->kt4,6)) 
      strncpy(phase_out,hdr->kt4,8);
    else
      strcpy(phase_out,"T4");
    return(hdr->t4); }
  else if (!strcasecmp(phase_name,"t5")) {
    if (strncasecmp("-12345",hdr->kt5,6)) 
      strncpy(phase_out,hdr->kt5,8);
    else
      strcpy(phase_out,"T5");
    return(hdr->t5); }
  else if (!strcasecmp(phase_name,"t6")) {
    if (strncasecmp("-12345",hdr->kt6,6)) 
      strncpy(phase_out,hdr->kt6,8);
    else
      strcpy(phase_out,"T6");
    return(hdr->t6); }
  else if (!strcasecmp(phase_name,"t7")) {
    if (strncasecmp("-12345",hdr->kt7,6)) 
      strncpy(phase_out,hdr->kt7,8);
    else
      strcpy(phase_out,"T7");
    return(hdr->t7); }
  else if (!strcasecmp(phase_name,"t8")) {
    if (strncasecmp("-12345",hdr->kt8,6)) 
      strncpy(phase_out,hdr->kt8,8);
    else
      strcpy(phase_out,"T8");
    return(hdr->t8); }
  else if (!strcasecmp(phase_name,"t9")) {
    if (strncasecmp("-12345",hdr->kt9,6)) 
      strncpy(phase_out,hdr->kt9,8);
    else
      strcpy(phase_out,"T9");
    return(hdr->t9); }
  else if (!strncasecmp(phase_name,hdr->kt0,8)) {
    return(hdr->t0); }
  else if (!strncasecmp(phase_name,hdr->kt1,8)) {
    return(hdr->t1); }
  else if (!strncasecmp(phase_name,hdr->kt2,8)) {
    return(hdr->t2); }
  else if (!strncasecmp(phase_name,hdr->kt3,8)) {
    return(hdr->t3); }
  else if (!strncasecmp(phase_name,hdr->kt4,8)) {
    return(hdr->t4); }
  else if (!strncasecmp(phase_name,hdr->kt5,8)) {
    return(hdr->t5); }
  else if (!strncasecmp(phase_name,hdr->kt6,8)) {
    return(hdr->t6); }
  else if (!strncasecmp(phase_name,hdr->kt7,8)) {
    return(hdr->t7); }
  else if (!strncasecmp(phase_name,hdr->kt8,8)) {
    return(hdr->t8); }
  else if (!strncasecmp(phase_name,hdr->kt9,8)) {
    return(hdr->t9); }
  else {
    warn_msg("Cannot find phase name match");
    return (float)NAN;
  }
}

void parse(int argc, char **argv, ms_params *par) {
  int iarg;
  char *dummy;
  logical read_data=FALSE, define_win=FALSE;
  /* Required arguments: initialise to 0 and check later */
  par->method=0;
  par->model=0;
  /* Default markers ( Defaults depend on other options ) */
  par->root[0]='\0';
  /* Defaults */
  strcpy(par->window.phase_start,"B");
  strcpy(par->window.phase_end,"E");
  par->window.offset_start=0.0;
  par->window.offset_end=0.0;
  par->window.taper=-1;
  par->make_grd=FALSE;
  par->dof_s=1;

  iarg=0;
  while(++iarg<argc) {
/*     VRB(printf("argv[%d]=%s\n",iarg,argv[iarg]);) */
    if(argv[iarg][0]!='-') {
      fprintf(stderr,"%s ",argv[iarg]);
      abort_msg("is not an option. This error can also be due to the wrong number of arguments for a previous option"); 
    }
    if(!strcasecmp(argv[iarg],"-data")) {
      read_data=TRUE;
      if ( iarg+2>=argc ) 
	abort_msg("-data option must be followed by two sac-filenames (horizontal components).");
      read_seis_file(argv[++iarg],&(par->hdr_hor1),&(par->data_hor1));
      read_seis_file(argv[++iarg],&(par->hdr_hor2),&(par->data_hor2));
 
      if (strlen(par->root)==0) {
	/* set default root name to root of first data file */
	if ( dummy=rindex(argv[iarg-1],'.') )
	  strncpy(par->root,argv[iarg-1],(size_t)(dummy-argv[iarg-1])); 
	else
	  strcpy(par->root,argv[iarg-1]);
      }
    }
    /* Method options */
    else if(!strcasecmp(argv[iarg],"-me") || !strcasecmp(argv[iarg],"-minevalue")) {
      par->method = MINEVALUE; }
    else if(!strcasecmp(argv[iarg],"-mt") || !strcasecmp(argv[iarg],"-mintransverse")) {
      par->method = MINTRANSVERSE; }
    else if(!strcasecmp(argv[iarg],"-conv")) {
      par->method = CONV; }
    else if(!strcasecmp(argv[iarg],"-cx") || !strcasecmp(argv[iarg],"-correl")) {
      par->method = CORREL; 
      if ( iarg+3>=argc ) 
	abort_msg("-correl option must be followed by maxshift (in seconds) and two sac-filenames (horizontal components).");
      par->method_q.cor_par.maxshift=atof(argv[++iarg]);
      read_seis_file(argv[++iarg],&(par->method_q.cor_par.hdr_ref1),&(par->method_q.cor_par.data_ref1));
      read_seis_file(argv[++iarg],&(par->method_q.cor_par.hdr_ref2),&(par->method_q.cor_par.data_ref2));
    }
    /* Splitting model options */
    else if(!strcasecmp(argv[iarg],"-single")) {
      par->model = SINGLE_HOR_SPLIT; 
      if ( iarg+3>=argc ) 
	abort_msg("-single option must be followed by 3 arguments (stepfast, stepdelay, maxdelay)");
      par->model_q.split_par.bot.faststep = atof(argv[++iarg]);
      par->model_q.split_par.bot.fastmin = 0;
      par->model_q.split_par.bot.fastmax = 180;
      par->model_q.split_par.bot.timestep = atof(argv[++iarg]);
      par->model_q.split_par.bot.timemin = 0;
      par->model_q.split_par.bot.timemax = atof(argv[++iarg]);

      par->model_q.split_par.top.faststep = par->model_q.split_par.top.timestep = 1.0;
      par->model_q.split_par.top.fastmin = par->model_q.split_par.top.timemin = 0.0;
      par->model_q.split_par.top.fastmax = par->model_q.split_par.top.timemax =0.0;
    }
    else if(!strcasecmp(argv[iarg],"-double")) {
      par->model = DOUBLE_HOR_SPLIT; 
      if ( iarg+3>=argc ) 
	abort_msg("-double option must be followed by 3 arguments (stepfast, stepdelay, maxdelay)");
      par->model_q.split_par.top.faststep=par->model_q.split_par.bot.faststep= atof(argv[++iarg]);
      par->model_q.split_par.top.fastmin =par->model_q.split_par.bot.fastmin = 0;
      par->model_q.split_par.top.fastmax =par->model_q.split_par.bot.fastmax = 180;
      par->model_q.split_par.top.timestep=par->model_q.split_par.bot.timestep= atof(argv[++iarg]);
      par->model_q.split_par.top.timemin =par->model_q.split_par.bot.timemin = 0;
      par->model_q.split_par.top.timemax =par->model_q.split_par.bot.timemax = atof(argv[++iarg]);
    }
    else if(!strcasecmp(argv[iarg],"-singlesub") ) {
      par->model = SINGLE_HOR_SPLIT; 
      if ( iarg+6>=argc ) 
	abort_msg("-singlesub option must be followed by 6 arguments (stepfast minfast maxfast stepdelay mindelay maxdelay)");
      par->model_q.split_par.bot.faststep = atof(argv[++iarg]);
      par->model_q.split_par.bot.fastmin = atof(argv[++iarg]);
      par->model_q.split_par.bot.fastmax = atof(argv[++iarg]);
      par->model_q.split_par.bot.timestep = atof(argv[++iarg]);
      par->model_q.split_par.bot.timemin = atof(argv[++iarg]);
      par->model_q.split_par.bot.timemax = atof(argv[++iarg]);

      par->model_q.split_par.top.faststep = par->model_q.split_par.top.timestep = 1.0;
      par->model_q.split_par.top.fastmin = par->model_q.split_par.top.timemin = 0.0;
      par->model_q.split_par.top.fastmax = par->model_q.split_par.top.timemax =0.0;
    }
    /* Window definition */
    else if(!strcasecmp(argv[iarg],"-winp") ) {
      if ( iarg+2>=argc ) 
	abort_msg("-winp must be followed by 2 arguments (two phase names)");
      define_win=TRUE;
      strncpy(par->window.phase_start,argv[++iarg],8); par->window.phase_start[8]='\0';
      par->window.offset_start=0.0;
      strncpy(par->window.phase_end,argv[++iarg],8); par->window.phase_end[8]='\0';
      par->window.offset_end=0.0;
    }
    else if(!strcasecmp(argv[iarg],"-wint") ) {
      if ( iarg+3>=argc ) 
	abort_msg("-winp must be followed by 3 arguments (phase name, offset_start,offset_end)");
      define_win=TRUE;
      strncpy(par->window.phase_start,argv[++iarg],8); par->window.phase_start[8]='\0';
      strcpy(par->window.phase_end,par->window.phase_start);
      par->window.offset_start=atof(argv[++iarg]);
      par->window.offset_end=atof(argv[++iarg]);
    }
    else if(!strcasecmp(argv[iarg],"-taper") ) {
      if ( iarg+1>=argc ) 
	abort_msg("-taper must be followed by 1 argument (taper length in sec)");
      par->window.taper=atof(argv[++iarg]);
    }
    else if(!strcasecmp(argv[iarg],"-grd") ) {
      if (!par->make_grd)
	par->make_grd=MAKE_GRD;
    }
    else if(!strcasecmp(argv[iarg],"-gmt") ) {
      par->make_grd=MAKE_GMT;
    }
    else if(!strcasecmp(argv[iarg],"-gmt5") ) {
      par->make_grd=MAKE_GMT5;
    }
    else if(!strcasecmp(argv[iarg],"-dof") ) {
      if ( iarg+1>=argc ) 
	abort_msg("-dof must be followed by 1 argument (degrees of freedom per sec)");
      par->dof_s=atof(argv[++iarg]);
    }
    else if(!strcasecmp(argv[iarg],"-name") ) {
      if ( iarg+1>=argc ) 
	abort_msg("-name must be followed by 1 argument (file name root)");
      strncpy(par->root,argv[++iarg],128);
      par->root[127]='\0';
    }
    else if(!strcasecmp(argv[iarg],"-v") ) {
      verbose=1;
    }
    else if(!strcasecmp(argv[iarg],"-h") ) {
      usage("multisplit");
    }
    else {
      fprintf(stderr,"%s ",argv[iarg]);
      abort_msg("is not a known option");
    }
  }
  /* Consistency and error checking */
  if (!read_data)
    abort_msg("Need to specify -data option");
  if(par->method==0)
    abort_msg("Need to choose a method (one of -minevalue -mintransverse -conv -correl)");
  if(par->model==0)
    abort_msg("Need to choose a splitting model (one of -single -double -singlesub)");
  /* set taper to 2 s if undefined and whole record is used */
  if (!define_win && par->window.taper<0 ) {
    warn_msg("Setting taper to 2 s for whole record mode.");
    par->window.taper=2.0;
  }    
}

void abort_msg(char *msg) {
  fprintf(stderr,"%s\n ABORT \n",msg);
  exit(10);
}

void warn_msg(char *msg) {
  fprintf(stderr,"WARNING: %s\n",msg);
}

void usage(char *cmd) {
  fprintf(stderr,"Usage: %s [OPTIONS] \n",cmd);
  fprintf(stderr,"\
OPTIONS:\n\
\n\
Required:\n\
\n\
-data hor1 hor2\n\
                         hor1 and hor2 are sac files containing the S phase for which splitting\n\
                         should be measured.  The header variables CMPAZ and BAZ must be set.\n\
\n\
Method: Choose one of\n\
\n\
-me or -minevalue        Minimise second eigenvalue (i.e. impose as linear a motion as possible)\n\
                         (Secondary result: polarisation) - e.g. Silver and Chan, 1991\n\
\n\
-mt or -mintransverse    Minimise transverse energy - e.g. Silver and Chan, 1991\n\
\n\
-conv                    Convolution method - Menke and Levin, 2003\n\
\n\
-cx or -correl maxshift refh1 refh2       Reference station method (see Eken and Tilmann BSSA, 2014, doi:10.1785/0120140020 )\n\
                         Minimise difference to reference station\n\
                         maxshift: maximum amount of travel time shift between stations\n\
                         refh1, refh2: sac file containing two horizontal components with\n\
                                       reference trace\n\
                         (Secondary results: time shift, scale factor)\n\
\n\
\n\
Splitting model: choose one of\n\
\n\
-single stepfast stepdelay maxdelay\n\
                         Single layer splitting. \n\
                         stepfast, stepdelay: stepwidths for grid search\n\
                         maxdelay: maximum splitting delay\n\
                         \n\
-double stepfast stepdelay maxdelay\n\
                         Two layer splitting\n\
\n\
-singlesub stepfast minfast maxfast stepdelay mindelay maxdelay \n\
                         Search selected fast directions\n\
\n\
Window definition:  Choose one of\n\
\n\
-winp t1 t2              n1,n2 can be arbitray header variable-names. If there are not one of the \n\
                         valid SAC header variables, multisplit looks for headernames matching those\n\
\n\
-wint t1 start end       Phase name or header variable and offset to start and end of analysis window\n\
\n\
DEFAULT                  Use the whole record.\n\
\n\
\n\
Optional modifiers\n\
\n\
-taper sec               Taper records with taper length sec.  If no tapering is carried out\n\
                         shift windows to calculate cross-correlation rather than calculating\n\
                         cross-correlation of windowed traces.\n\
\n\
-grd                     Convert 2-parameter error surface binary files to GMT grid files\n\
                                        \n\
-gmt                     Create and execute GMT script for display - write out time series files (GMT4)\n\
-gmt5                     Create and execute GMT script for display - write out time series files (GMT5)\n\
\n\
-name root               Set root of output file names\n\
                         (Default: derive name from first input file root)\n\
\n\
-dof f                   Set number of degrees of freedom per second (Default: 1s)\n\
\n\
-v                       Verbose output\n\
\n\
Other options\n\
\n\
-h                       Show this help text\n\
\n\
");
  exit(10);
}

  /* Spares:
                         The first digit gives the parameter number for x-axis, the second 
                         digit the parameter for the y-axis (e.g. -grd 34 for the two
                         layer problem would give the bottom layer parameters)
  */
