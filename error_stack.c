/* errror_stack */
/* Author: F Tilmann */
/* Contact: tilmann|a|gfz-potsdam.de */

/* Stack error surface files produced by multisplit and find best overall
   splitting parameter */


/* (C) 2004 F Tilmann */
/* This code is released under the GNU public license */

/* Code uses gsl and gslblas library */
/* History: */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stddef.h>
#include <math.h>
#include <errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


int verbose=0;
#define VRB(command) { if(verbose) { command  ; fflush(stdout); }}
#define ASSERT(cond,msg) { if (!(cond)) {   fprintf(stderr,"ASSERTION VIOLATION: %s\n ABORT \n",msg);  exit(10);}}
#define WRITEVEC(name,vec) { {\
	FILE *vectorout;\
	vectorout=fopen(name,"w");\
	if (!vectorout) { \
      sprintf(abort_str,"Cannot open %s for output.",name);\
      abort_msg(abort_str);\
    }\
    gsl_vector_fprintf(vectorout,vec,"%f");\
    fclose(vectorout); }\
  }

#if !defined NAN
#define NAN (strtod("NAN",NULL))
#endif

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))
/* PER(per,x) return a value v with range 0..per-1 such that v+n*per=x (n an integer). It only works for integers */
#define PER(per,x) ( x<0 ? (x + per*(1+(-x)/per)) : (x%per) )


#define MAXFILES 256
#define MAXDIM 4



typedef struct {
  int gmt;
  int weight;
  char root[256];
  int nfiles;
  int mcsamples;
  int bootstrap_samples; // Number of bootstrap samples
  int use_exact; 
  float scale_dof;
  char *fnames[MAXFILES];
} params;

const gsl_rng_type * rng_type;
gsl_rng * rng;


char warn_str[1024];
char abort_str[1024];

FILE *open_for_write(char *root,char *extension);
FILE *open_for_read(char *root,char *extension);
void parse(int argc, char **argv, params *par);
void warn_msg(char *msg);
void abort_msg(char *msg);
void usage(char *cmd);
double invfisher(double nu1,double nu2,double conf);
double betai(double a, double b, double x);


void ind2sub(int *sub,long int ind1d, int sizes[], int rank) {
  /* convert one-dimensional index to multi-dimensional indices */
  /* Input:   int1d  : 1D index
     sizes[] : Arrays giving the number of elements in each dimension
     rank: number of dimensions (indices)
     Output:  sub[]   A vector containing the multi-dimensional indices */
  int i;

  for(i=rank-1; i>=0; i--) {
    sub[i]=ind1d % sizes[i];
    /*       imin_1d-=imin[k][i]; */ // redundant as this part is removed by integer division in next line anywat 
    ind1d /= sizes[i]; 
    // VRB(printf("DEBUG i: %d imin_1d: %ld m[i]: %d  imin[i]: %d\n",i,imin_1d,m[i],imin[k][i]))
  }
  ASSERT(ind1d==0,"ind1d should be zero at this stage, as all indices determined and there is no larger dimension left.")
}

long *bootstrap_analysis(gsl_vector *v_err_surf[],int nfiles,int bootstrap_samples, params *par) {
  int l,k; // loop variables 
  unsigned long int r;
  long imin;
  int use_count[nfiles]; // how many times is a file represented in a bootstrap sample
  float weight=1;  // as default use 1 for weight
  long *imin1d_bootstrap=malloc(bootstrap_samples*sizeof(long)); 

  gsl_vector *v_err_stack=gsl_vector_alloc(v_err_surf[0]->size);
  gsl_vector *v_temp=gsl_vector_alloc(v_err_surf[0]->size);

  for (l=0; l<bootstrap_samples; l++) { // Loop over bootstrap samples
    VRB(printf("\rbootstrap sample %6d",l));
    // reset use_counter 
    for(k=0;k<nfiles;k++) 
      use_count[k]=0;
    // reset stacked error srufcae
    gsl_vector_set_zero(v_err_stack);


    // generate bootstrap sample counts for each file
    for(k=0;k<nfiles;k++) {
      r=gsl_rng_uniform_int(rng, nfiles);
      use_count[r]++;
    }

    for(k=0;k<nfiles;k++) { /* apply multipliers */
      if (use_count[k]==0) 
	continue;          // can skip rest if event is not used in this bootstrap samples
      if (par->weight) {
	// weight with the inverse of minimum energy (if relative weighting is selected)
	imin=gsl_vector_min_index(v_err_surf[k]);
	weight=1/gsl_vector_get(v_err_surf[k],imin);
      } 
      gsl_vector_memcpy(v_temp,v_err_surf[k]);
      gsl_vector_scale(v_temp,weight);
      gsl_vector_add(v_err_stack,v_temp);
    }
    // find index of minimum in bootstrapped stacks
    imin1d_bootstrap[l]=gsl_vector_min_index(v_err_stack);
  }
  // free allocated vector space
  gsl_vector_free(v_err_stack);
  gsl_vector_free(v_temp);
  return(imin1d_bootstrap);
}



int main(int argc, char **argv)
{
  params *par=(params *) malloc(sizeof(params));
  FILE *hdr_file, *bin_file;
  FILE *output;
  int i,j,k;
  /* Properties of each file */
  float split_par;
  char methodstring[128],label[MAXDIM][128];
  int dim,m[MAXDIM];
  float min[MAXDIM],max[MAXDIM],step[MAXDIM];
/*   gsl_matrix *err_surf; */
  gsl_matrix *err_stack;  // for 2D error-surfaces

  gsl_vector *v_err_surf[MAXFILES], *v_err_stack, *v_temp;  // for arbitrary dimension error surfaces
  float weight;
  /* Properties of first file */
  char rmethodstring[128],rlabel[MAXDIM][128];
  int rdim,rm[MAXDIM];
  float rmin[MAXDIM],rmax[MAXDIM],rstep[MAXDIM];
  /* Variables for ensemble (including arrays for remembering a value for each file) */
  int imin[MAXFILES][MAXDIM],imint[MAXDIM],imcmc[MAXDIM],ibootstrap[MAXDIM];
  long imin_1d,imint_1d,tot_length;
  float dof[MAXFILES],postconf[MAXFILES];
  float tot_dof,valmin[MAXFILES];
  float tot_weight;
  long *imin1d_bootstrap;

  double maxval,emin,misfit,logprob,value,best[MAXDIM],lbound[MAXDIM],ubound[MAXDIM],err[MAXDIM];
  /* error analysis */
  double conf[9]={.68,  .95,.99, .999,.9999,.99999,.999999,.9999999, .99999999 };  /* conf[1] is the level of significance */
  double contour[9];
  double null;
  int j1,j2,k1,k2;
  int periodic;
  char cmdstring[1024];
  gsl_vector_view dvue1;
  gsl_matrix_view vue_matrix;

  int status;
  char rejectstring[128];

  // set up global variables related to random number generation
  gsl_rng_env_setup();
  rng_type = gsl_rng_default;
  rng = gsl_rng_alloc (rng_type);

  // parse command line arguments
  parse(argc,argv, par);

  tot_dof=0;
/*   tot_weight=0; */

  for (k=0;k<par->nfiles;++k) {  /* read loop */
    hdr_file=open_for_read(par->fnames[k],".hdr");
    bin_file=open_for_read(par->fnames[k],".bin");
    if ( fscanf(hdr_file,"%s %f %f\n",methodstring,&split_par,&dof[k]) != 3 ) {
      fprintf(stderr,"Format error in %s.hdr file line 1",par->fnames[k]); 
      abort_msg(abort_str);
    }
    if ( fscanf(hdr_file,"%d\n",&dim) != 1 ) {
      sprintf(abort_str,"Format error in %s.hdr file line 2",par->fnames[k]); 
      abort_msg(abort_str);
    }
    if (dim>MAXDIM) 
      abort_msg("Maximum number of dimensions exceeded. Increase MAXDIM in source and recompile");
    for (i=0; i<dim;i++) {
      if ( fscanf(hdr_file,"%d",&m[i]) != 1 ) {
	sprintf(abort_str,"Format error in %s.hdr file line 3",par->fnames[k]); 
	abort_msg(abort_str);
      }
    }
    for (i=0; i<dim;i++) {
      if (fscanf(hdr_file,"%s %f %f %f\n",label[i],&min[i],&max[i],&step[i]) != 4) {
	sprintf(abort_str,"Format error in %s.hdr file line %d",par->fnames[k],4+i); 
	abort_msg(abort_str);
      }
    }
    if(k==0) {
      strcpy(rmethodstring,methodstring);
      rdim=dim;
      tot_length=1;
      for(i=0;i<dim;i++) {
	rm[i]=m[i]; rmin[i]=min[i]; rmax[i]=max[i]; rstep[i]=step[i];
	strcpy(rlabel[i],label[i]);
	tot_length*=m[i];
      }
      v_err_stack=gsl_vector_calloc(tot_length);
      v_temp=gsl_vector_alloc(tot_length);
/*       if (rdim==2) { */
/* 	err_stack=gsl_matrix_alloc(rm[0],rm[1]); */
/* 	gsl_matrix_set_zero(err_stack);  */
/* 	err_surf =gsl_matrix_alloc(rm[0],rm[1]); */
/*       } else { */
/* 	abort_msg("Currently code not set up for number of dimensions not equal 2"); */
/*       }			    */
    } else { /* k!=0: check that definitions consistent with first file */
      if (strcmp(rmethodstring,"Mixed") && strcmp(rmethodstring,methodstring)) {
	sprintf(warn_str,"Method in %s and %s disagree:\n %s vs %s\n This is OK if you wish to combine different methods",par->fnames[0],par->fnames[k],rmethodstring,methodstring);
	strcpy(rmethodstring,"Mixed");
	warn_msg(warn_str);
      }
      if (rdim != dim) {
	sprintf(abort_str,"Number of dimensions in %s and %s disagree:\n %d vs %d",par->fnames[0],par->fnames[k],rdim,dim);
	abort_msg(abort_str);
      }
      for(i=0;i<dim;i++) {
	if(rm[i]!=m[i] || rmin[i]!=min[i] || rmax[i]!=max[i] || rstep[i]!=step[i] || strcmp(rlabel[i],label[i])) {
	sprintf(abort_str,"There is a discrepancy between %s and %s disagree in the definition of dimension %d\n",par->fnames[0],par->fnames[k],i);
	abort_msg(abort_str);
	}
      }
    }
    // dimension agnostic code 
    v_err_surf[k]=gsl_vector_alloc(tot_length);
    if( gsl_vector_fread(bin_file,v_err_surf[k]) ) {
      sprintf(abort_str,"There was a problem reading %s.bin",par->fnames[k]);
      abort_msg(abort_str);
    }
  } /* end of read loop */

  for (k=0;k<par->nfiles;++k) { /* Analysis loop */
    imin_1d=gsl_vector_min_index(v_err_surf[k]);
    valmin[k]=gsl_vector_get(v_err_surf[k],imin_1d);
/*    Code specific to 2D matrix
/*     if( gsl_matrix_fread(bin_file,err_surf) ) { */
/*       sprintf(abort_str,"There was a problem reading %s.bin",par->fnames[k]); */
/*       abort_msg(abort_str); */
/*     } */
//     /* find minimum index  */
/*     gsl_matrix_min_index(err_surf,(size_t *)&imin[k][0],(size_t *)&imin[k][1]); */
/*     valmin[k]=gsl_matrix_get(err_surf,imin[k][0],imin[k][1]); */
/*     printf("%25s %s: %6f  %s: %6f emin: %f dof: %f\n",par->fnames[k], */
/* 	   rlabel[0],min[0]+imin[k][0]*step[0], */
/* 	   rlabel[1],min[1]+imin[k][1]*step[1],valmin[k],dof[k]); */
    dof[k]*=par->scale_dof;
    printf("%25s emin:%f  dof: %f",par->fnames[k],valmin[k],dof[k]);
    //    VRB(printf("\n DEBUG imin_1d: %ld  dim: %d\n",imin_1d,dim));
    ind2sub(&imin[k][0],imin_1d,rm,dim);

    for(i=0; i<dim; i++) {
      printf(" %s: %6f", rlabel[i],min[i]+imin[k][i]*step[i]);
    }
    printf("\n");
    if (par->weight)  /* normalise by minimum value if requested */
      weight=1/valmin[k];
    else
      weight=1;
/* Code for 2D matrix */
/*     gsl_matrix_scale(err_surf,weight); */
/*     gsl_matrix_add(err_stack,err_surf); */
    /* dimension agnostic */
    gsl_vector_memcpy(v_temp,v_err_surf[k]);
    gsl_vector_scale(v_temp,weight);
    gsl_vector_add(v_err_stack,v_temp);


    tot_dof += weight*(dof[k]+par->scale_dof*split_par);;    
      /* I am not sure if this way of manipulating degrees of reedom is correct.*/
      /* However, we must introduce some weighing of the degrees of freedom, otherwise */
      /* I could take one good measurement and lots of bad ones, get the results of the */
      /* good one but with apparently much diminished error because of the greater number */
      /* of DOFs */
     tot_weight+= weight; 
  } /* end of analysis loop */ 


  /* Analysis of error surface */
  /*  2D code: */
  /*   gsl_matrix_min_index(err_stack,(size_t *)&imint[0],(size_t *)&imint[1]); */
  /*   emin=gsl_matrix_get(err_stack,imint[0],imint[1]); */
  imint_1d=gsl_vector_min_index(v_err_stack);
  emin=gsl_vector_get(v_err_stack,imint_1d);
  ind2sub(imint,imint_1d,rm,dim);

  for (i=0;i<dim;i++) {
    best[i]=rmin[i]+imint[i]*step[i];
  }

  tot_dof=par->nfiles*tot_dof/tot_weight;
  tot_dof-=par->scale_dof*split_par;
  

  VRB(printf("DEBUG: before calculating confidence intervals\n"));

  float tot_dof_eff=tot_dof;
  float tot_dof_thresh=200.;
  if ( tot_dof>tot_dof_thresh) {
    fprintf(stderr,"WARNING: for very large total number of degrees of freedom \n\
(equiv. total length of analysed data), the confidence interval calculation leads \n\
to underflow. The effective number of degrees of freedom for conf. level calculation\n\
 has thus been reduced from %.2f to %.0f\n",tot_dof,tot_dof_thresh);
    tot_dof_eff=tot_dof_thresh;
  } 
  for(i=0;i<9;i++) { /* loop over confidence intervals */
    VRB(printf("i: %d confidence %f  tot_dof %f split_par %f\n",i,conf[i],tot_dof,split_par));
    contour[i]=1+split_par*invfisher((double)split_par,(double)tot_dof_eff,conf[i])/tot_dof_eff;
    VRB(printf("i: %d confidence %f contour %f tot_dof %f split_par %f\n",i,conf[i],contour[i],tot_dof,split_par));
  }
  VRB(printf("DEBUG: after calculating confidence intervals\n"));

  printf("Number of files in stack:   %d\n",par->nfiles);
  printf("Total degrees of freedom:   %f\n",tot_dof);
  printf("Norm. energy minimum:       %f\n",emin/tot_weight);

  if ( dim==2 ) {
    vue_matrix=gsl_matrix_view_vector(v_err_stack,m[0],m[1]);
    err_stack=&vue_matrix.matrix;
    if (par->gmt) {
      output=open_for_write(par->root,".cont");
      for (i=0;i<9;i++){
	fprintf(output,"%f %s\n",contour[i]*emin, (i==1 ? "A" : "C" ));
      }
      fclose(output);
    }
    null=-1;
    for(i=0;i<9;i++) {
      if (rmin[1]==0.0 && contour[i]*emin<gsl_matrix_get(err_stack,0,0))
	    null=conf[i];
    }

    if ((null)<0) 
      strcpy(rejectstring, "UNK");
    else if (null==0.0)
      strcpy(rejectstring, "<68");
    else
      sprintf(rejectstring,"%7.4f",null*100);
    fprintf(stderr,"DEBUG: Null %f Rejectstring |%s|\n",null,rejectstring);
  
    /* Global search parameter 1 (normally this is splitting delay): */
    for (k1=rm[1]-1; k1>=imint[1];k1--) {
      dvue1=gsl_matrix_column(err_stack,k1);
      if (gsl_vector_min(&dvue1.vector) <= contour[1]*emin) 
	break;
    }
    for (k2=0; k2<=imint[1];k2++) {
      dvue1=gsl_matrix_column(err_stack,k2);
      if (gsl_vector_min(&dvue1.vector) <= contour[1]*emin) 
	break;
    }

    lbound[1]=rmin[1]+k2*rstep[1];
    ubound[1]=rmin[1]+k1*rstep[1];
    err[1]=MAX(ubound[1]-best[1],best[1]-lbound[1]);
    VRB(printf("Global search time %f-%f (E: %f)\n",lbound[1],ubound[1],err[1]));
    if (k1==rm[1]-1 && k2==0)     /* both bounds at limit of grid search -> no constraints */
      /*     err_time=nan("");    */  /* set to NaN */ 
      err[1]=NAN;
    else if (k1==rm[1]-1 )     /* upper bounds at limit of grid search */
      err[1]=-err[1];    /* set error to negative (as flag) */

    /* Global search parameter 0: (normally this is the fast direction */
    if(!strcasecmp(rlabel[0],"Fast") && 180+rmin[0]-rmax[0]<2*rstep[0]) {
      VRB(printf("Periodic fast direction (parameter 0)\n"));
      periodic=1;
    } else {
      VRB(printf("Non-Periodic fast direction (parameter 1)\n"));
      periodic=0; 
    }
    for (j1=(periodic ? imint[0]+rm[0]/2 : rm[0]-1); j1>=imint[0];j1--) {
      dvue1=gsl_matrix_row(err_stack,PER(rm[0],j1));
      if (gsl_vector_min(&dvue1.vector) <= contour[1]*emin) 
	break;
    }
    for (j2=(periodic ? imint[0]-rm[0]/2 : 0); j2<=imint[0];j2++) {
      dvue1=gsl_matrix_row(err_stack,PER(rm[0],j2));
      if (gsl_vector_min(&dvue1.vector) <= contour[1]*emin) 
	break;
    }
    ubound[0]=rmin[0]+j1*step[0];
    lbound[0]=rmin[0]+j2*step[0];
    err[0]=MAX(ubound[0]-best[0],best[0]-lbound[0]);
    VRB(printf("Global search fast direction %f-%f (E: %f)\n",lbound[0],ubound[0],err[0]));

    if (err[0]>75. && periodic)
      err[0]=NAN;
    else if (!periodic && j1==rm[0]-1 && j2==0)     /* both bounds at limit of grid search -> no constraints */
      /*     err_time=nan("");    */  /* set to NaN */ 
      err[0]=NAN;
    else if (!periodic && ( j1==rm[0]-1 || j2==0 )  )     /* upper bounds at limit of grid search */
      err[0]=-err[1];    /* set error to negative (as flag) */
    printf("%24s:    %3.0f  +- %3.0f  (% 4.0f - %3.0f )\n",rlabel[0],best[0],err[0],lbound[0],ubound[0]);
    printf("%24s:    %4.2f +- %4.2f (% 4.2f - %4.2f )\n",rlabel[1],best[1],err[1],lbound[1],ubound[1]); 
    if (strncmp(rejectstring,"UNK",3))
      printf("Reject Null (%%) : %s (E: %f)\n",rejectstring,gsl_matrix_get(err_stack,0,0)/tot_weight);
  }   /* end of special part for single layer splitting  */
  else {
    /* part for dim>2, no more plotting, no more upper, lower bound, instead draw from probability distribution */
/*     gsl_vector *v_logprob_surf =v_err_surf;   // resuse already allocated vector to save memory but give it a new name for transparency */
    gsl_vector *v_logprob_surf =gsl_vector_calloc(tot_length);
    // display best solution
    for (i=0;i<rdim;i++) 
      printf("%24s: %f \n",rlabel[i],best[i]);
    // convert everything into logprobabilities

    for (i=0;i<v_logprob_surf->size; i++) {
	misfit=gsl_vector_get(v_err_stack,i);
        // use distribution for unknown error level and number of data points equal to number of Degrees of freedom       
        // Presumably because the number of data points is very large, this puts nearly all probability onto the minimum misfit 
        // point 
	logprob=-((tot_dof+split_par)/2) * log(misfit);   // pow(misfit,-(tot_dof+split_par)/2.);
        // Alternative: assume error level given by minimum
        // logprob=exp(-misfit/(2*emin));
	gsl_vector_set(v_logprob_surf,i,logprob); 
    }

    // Convert to absolute probabilities; we need to substract maximum of logarithm to avoid overflow
    maxval=gsl_vector_max(v_logprob_surf);
    for (i=0;i<v_logprob_surf->size; i++) { 
        value=exp(gsl_vector_get(v_logprob_surf,i)-maxval);
	gsl_vector_set(v_logprob_surf,i,value);
    }
    // WRITEVEC("probabilities.txt",v_logprob_surf);  
    // WRITEVEC("misfit.txt",v_err_stack);  
    
    // Calculate cumulative sum 
    for (i=1;i<v_logprob_surf->size; i++) {
	value=gsl_vector_get(v_logprob_surf,i-1)+gsl_vector_get(v_logprob_surf,i);
	gsl_vector_set(v_logprob_surf,i,value);
    }
    // WRITEVEC("cumprob.xy",v_logprob_surf);  
    
    // generate Monte-Carlo samples
    maxval=gsl_vector_get(v_logprob_surf,i-1);
    VRB(printf("Generating MC samples\n"));

    output=open_for_write(par->root,"_mc.x");
    fprintf(output,"#");
    for (j=0;j<rdim;j++) 
      fprintf(output,"%s ",rlabel[j]);
    fprintf(output,"\n");
    for (i=0;i<par->mcsamples;i++) {
        double r;
	r=gsl_rng_uniform(rng)*maxval;
        for (k=0;k<v_logprob_surf->size; k++) {
	    if (r<gsl_vector_get(v_logprob_surf,k)) 
	       break;
	}
	ind2sub(imcmc,k,rm,dim);
	for (j=0;j<dim;j++) {
            double out=rmin[j]+imcmc[j]*step[j];
            // randomize within each box covered
            if (!par->use_exact) {
	       out+= (-0.5+gsl_rng_uniform(rng))*step[j];
	    }
	    fprintf(output,"%f ",out);
	}
	fprintf(output,"\n");
    }
    fclose(output);
  }

  if (par->bootstrap_samples>0) {
    VRB(printf("Generating bootstrap samples\n"));
    imin1d_bootstrap=bootstrap_analysis(v_err_surf,par->nfiles,par->bootstrap_samples,par);
    // output bootstrap samples

    output=open_for_write(par->root,"_bootstrap.x");
    fprintf(output,"#");
    for (j=0;j<rdim;j++) 
      fprintf(output,"%s ",rlabel[j]);
    fprintf(output,"\n");
    for (i=0;i<par->bootstrap_samples;i++) {
 	ind2sub(ibootstrap,imin1d_bootstrap[i],rm,dim);
	for (j=0;j<dim;j++) {
            double out=rmin[j]+ibootstrap[j]*step[j];
            // randomize within each box covered
            if (!par->use_exact) {
	       out+= (-0.5+gsl_rng_uniform(rng))*step[j];
	    }
	    fprintf(output,"%f ",out);
	}
	fprintf(output,"\n");
    }
    fclose(output);
  }

  /* write bin file */
  output=open_for_write(par->root,".hdr");
  fprintf(output,"Stack_%s %f %f\n",rmethodstring,split_par,tot_dof);
  fprintf(output,"%d\n",rdim);
  for(i=0;i<rdim;i++)
    fprintf(output,"%d ",rm[i]);
  fprintf(output,"\n");
  for(i=0;i<rdim;i++)
    fprintf(output,"%-8s %f %f %f\n",rlabel[i],rmin[i],rmax[i],rstep[i]);
  fclose(output);
  output=open_for_write(par->root,".bin");
/*   gsl_matrix_fwrite(output,err_stack); */
  gsl_vector_fwrite(output,v_err_stack);
  fclose(output);

  /* now go through all input error surfaces again, and check which confidence interval the best choice corresponds to */
  for (k=0;k<par->nfiles;++k) {
    double x;
    // no need to re-read as all inputs are stored in memory
/*     bin_file=open_for_read(par->fnames[k],".bin"); */
//     gsl_matrix_fread(bin_file,err_surf);  /* we read this successfully before */ 
//     emin=gsl_matrix_get(err_surf,imint[0],imint[1]); /* get energy value at best splitting parameters of ensemble */ 
/*     gsl_vector_fread(bin_file,v_err_surf); */
    emin=gsl_vector_get(v_err_surf[k],imint_1d);
    VRB(printf("File %d: %s emin %f  valmin[k]: %f\n",k,par->fnames[k],emin,valmin[k]));

    /* this is the forward calculation of the confidence level from the energy ratio */
    x=(emin/valmin[k]-1)*dof[k]/split_par;
    x=MIN(1,(double)(dof[k]/(dof[k]+split_par*x))); /* for exact co-incidence round-off can lead to x values slightly larger than 1, hence have to use MIN */
    VRB(printf("X=%f\n",x));
    postconf[k]=1-betai((double)dof[k]/2.,(double)split_par/2,x);
    printf("Confidence value for %s : %f\n",par->fnames[k],postconf[k]);
/*     fclose(bin_file); */
  }
  
  if (par->gmt && dim==2) {
    /*convert bin error surface to grd file*/
    sprintf(cmdstring,"xyz2grd %s.bin -D%s/%s/%s/1/0/SingleSplit/\"Created by error_stack\" -G%s.grd -I%f/%f -R%f/%f/%f/%f -ZBLd",
	    par->root,rlabel[1],rlabel[0],methodstring,par->root,
	    rstep[1],rstep[0],
	    rmin[1],rmax[1],rmin[0],rmax[0]);
    VRB(printf("GRD conversion: %s\n",cmdstring));
    status=system(cmdstring);
    if (status) {
      sprintf(warn_str,"External GMT command xyz2grd execution failed. Error status: %d",status);
      warn_msg(warn_str);
    }
    output=open_for_write(par->root,".gmt");
    fprintf(output,"#!/bin/csh\n");
    fprintf(output,"# script auto-generated by multisplit\n");
    fprintf(output,"\n");
    fprintf(output,"# Variable parameters:\n");
    fprintf(output,"set root=%s\n",par->root);
    fprintf(output,"set best0=%f\n",best[0]);
    fprintf(output,"set best1=%f\n",best[1]);
    fprintf(output,"set label0=\"%s\"\n",rlabel[0]);
    fprintf(output,"set label1=\"%s\"\n",rlabel[1]);
    fprintf(output,"cat > /tmp/$root.description <<EOF\n");
    fprintf(output,"> 10 29 14 0 0 CT 0.564 20  c\n");
    fprintf(output,"%s %s stack of %d files.\n\n",rmethodstring,(par->weight ? "weighted" : "unweighted" ),par->nfiles);
    fprintf(output,"Energy Minimum %f (normalised: %f) DOF %f\n\n",emin,emin/tot_weight,tot_dof);
    fprintf(output," %s %3.0f \\261 %3.0f %s %4.2f \\261 %4.2f RejectNull %s\n",
	    rlabel[0],best[0],err[0], rlabel[1],best[1],err[1],rejectstring);
    fprintf(output,"EOF\n\n");
    fprintf(output,"cat > /tmp/${root}_allsplit.xy <<EOF\n");
    for (k=0;k<par->nfiles;k++) 
      fprintf(output,"%f %f\n",min[1]+imin[k][1]*step[1],min[0]+imin[k][0]*step[0]);
    fprintf(output,"EOF\n\n");
    fprintf(output,"\
### Everything below this line is independent of the particular event used\n\
 \n\
set grdrange=`grdinfo -C ${root}.grd | awk '{print $2 \"/\" $3 \"/\" $4 \"/\" $5 }'`\n\
set psfile=${root}.ps\n\
 \n\
gmtdefaults -D >.gmtdefaults\n\
gmtset PAGE_ORIENTATION portrait MEASURE_UNIT cm WANT_EURO_FONT TRUE LABEL_FONT_SIZE 12 ANOT_FONT_SIZE 10 PAPER_MEDIA a4 D_FORMAT %%lg\n\
\n\
# 3cm Descriptive text\n\
pstext -M -X0 -Y0 -R0/20/0/29 -Jx1 -K > $psfile </tmp/${root}.description\n\
 \n\
# 8 cm Error surface\n\
grdcontour -X2 -Y20.5 ${root}.grd -C${root}.cont -R$grdrange -JX17/6.5 -B0.5:\"$label1\":/20:\"$label0\":WSen -O -K -A-1f1 -G1000 -Wa1.5p -Wc0.5p >>$psfile\n\
psxy -R -JX -Sx0.3 -W2p/200/100/100  -O -K  >>$psfile <<EOF\n\
$best1 $best0\n\
EOF\n\
psxy < /tmp/${root}_allsplit.xy -R -JX -Sx0.2 -W1p  -O -K  >>$psfile \n\
\n\
psxy < /dev/null -Jx1 -R -O >>$psfile\n\
\\rm /tmp/${root}_allsplit.xy /tmp/${root}.description\n\
");
    fclose(output);
    sprintf(cmdstring,"csh %s.gmt",par->root);
    VRB(printf("Executing GMT script %s\n",cmdstring));
    status=system(cmdstring);
/*     status=system("echo hello world; csh err_stack_test.gmt ; echo hello world again"); */
    if (status) {
      sprintf(warn_str,"Execution of GMT script failed. Error status: %d",status);
      warn_msg(warn_str);
    }
  }
  return 0;
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

FILE *open_for_read(char *root,char *extension){
  char tmpstring[256];
  FILE *output;
  strcpy(tmpstring, root);
  strcat(tmpstring, extension);
  
  output=fopen(tmpstring,"r");
  if (!output) { 
    sprintf(abort_str,"Cannot open %s for input.",tmpstring);
    abort_msg(abort_str);
  }
  return(output);
}


void parse(int argc, char **argv, params *par) {
  int iarg;
  char *dummy;
  /* Default markers ( Defaults depend on other options ) */
  par->root[0]='\0';
  par->weight=0;
  par->gmt=0;
  par->nfiles=0;
  par->use_exact=0;
  par->mcsamples=100;
  par->bootstrap_samples=0; 
  par->scale_dof=1.;

  iarg=0;
  while(++iarg<argc) {
    if(argv[iarg][0]!='-') {
      /* argument not an option */
      break; 
    }
    if(!strcasecmp(argv[iarg],"--")) {
      iarg++;
      break; } 
    else if(!strncasecmp(argv[iarg],"-exact",6)) {
      par->use_exact = 1; }
    else if(!strncasecmp(argv[iarg],"-weig",5)) {
      par->weight = 1; }
    else if(!strncasecmp(argv[iarg],"-gmt",4)) {
      par->gmt = 1; }
    else if(!strncasecmp(argv[iarg],"-name",5) ) {
      if ( iarg+1>=argc ) 
	abort_msg("-name must be followed by 1 argument (file name root)");
      strncpy(par->root,argv[++iarg],256); par->root[255]='\0'; 
    }
    else if(!strncasecmp(argv[iarg],"-mc",3) ) {
      if ( iarg+1>=argc ) 
	abort_msg("-mc must be followed by 1 argument (number of samples)");
      par->mcsamples=(int)strtol(argv[++iarg],NULL,10); 
      if (errno || par->mcsamples<=0 ) {
	 abort_msg("Argument of -mc must be a positive integer number");
      }
    }
    else if(!strncasecmp(argv[iarg],"-bootstrap",10) ) {
      if ( iarg+1>=argc ) 
	abort_msg("-bootstrap must be followed by 1 argument (number of samples)");
      par->bootstrap_samples=(int)strtol(argv[++iarg],NULL,10); 
      if (errno || par->bootstrap_samples<0 ) {   // accept 0 as meaning effectively no bootstrap samples
	 abort_msg("Argument of -bootstrap must be a positive integer number");
      }
    }
    else if(!strncasecmp(argv[iarg],"-scale-dof",10) ) {
      if ( iarg+1>=argc )  
	abort_msg("-scale-dof must be followed by 1 argument (float)");
      par->scale_dof=(float)strtod(argv[++iarg],NULL);
      if (errno || par->scale_dof<=0 ) {
	 abort_msg("Argument of -mc must be a positive (float) number");
      }
    }
    else if(!strncasecmp(argv[iarg],"-h",2)) {
      usage("error_stack"); }
    else if(!strncasecmp(argv[iarg],"-v",2)) {
      verbose=1; }
    else {
      fprintf(stderr,"%s ",argv[iarg]);
      abort_msg("is not a known option");
    }
  }
  for (iarg=iarg;iarg<argc;++iarg) {
    par->fnames[par->nfiles]=strdup(argv[iarg]);
    dummy=rindex(argv[iarg],'.');
    if (!strcmp(dummy,".bin")) {
      /* strip .bin extension */
      par->fnames[par->nfiles][(size_t)(dummy-argv[iarg])]='\0';
    }
    if (strlen(par->root)==0) {
      strncpy(par->root,par->fnames[par->nfiles],240); 
      strcat(par->root,"_stack");
    }
    ++par->nfiles;
    if (par->nfiles >= MAXFILES ) {
      sprintf(abort_str,"ERROR: maximum number of files (%d) exceeded. Increase MAXFILES in error_stack.c and recompile\n",MAXFILES);
      abort_msg(abort_str);
    }
  }
  if (par->nfiles==0)
    abort_msg("For stacking we need at least one surface to stack.  Use -h option to get help");
}

void abort_msg(char *msg) {
  fprintf(stderr,"%s\n ABORT \n",msg);
  exit(10);
}

void warn_msg(char *msg) {
  fprintf(stderr,"WARNING: %s\n",msg);
}

void usage(char *cmd) {
  fprintf(stderr,"Usage: %s [OPTIONS] file1 file2 file3 ... \n",cmd);
  fprintf(stderr,"\
file1 ... must be .bin files (output of e.g. multisplit). \n\
The corresponding .hdr files must also be present.  Note that the .bin\n\
extension of the files can be omitted.\n\
\n\
OUTPUT FILES\n\
\n\
.bin, .hdr            the stacked error surfaces\n\
\n\
.grd,.gmt,.ps         grdfile, GMT scripts and postscript for visualisation \n\
                      of result (only for number of parameters equal 2)\n		\
                       (if -gmt option has been selected)\n\
\n\
OPTIONS:\n\
\n\
-weight               normalise error surfaces by minimum value before\n\
                      stacking (default is no normalisation)\n\
\n\
-name <root>           Set root of output file names\n\
                      (Default: derive name from first input file root)\n\
\n\
-gmt                  Plot results with GMT (only for number of parameters equal 2\n	\
\n\
-mc <samples>         How many Monte Carlo samples to use (only for number of parameters not equal two\n\
                      [ Default 100]\n\
-exact                When calculating sample of distribution, use the exact point of the grid search\n\
                      [Default: randomize sample point within half a step length either way of the gridpoint]\n\
-scale-dof <scalefactor> Used to scale the input degrees-of-freedom. Overestimate in number of degrees of\n\
                      results in underestimate in uncertainty. So if the uncertainty appears too small,\n\
                      use a scalefactor less than 1, if they appear too big, use a scalefactor of more than 1.\n\
\n\
-bootstrap <samples>   Use bootstrap over error surfaces to estimate uncertainty\n\
                      Also affected by -exact option\n\
                      [ Default do not estimate bootstrap errors]\n\
\n\
-v                    Verbose output\n\
");
  exit(10);
}
