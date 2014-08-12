#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector_float.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include "gsl_seis.h"

#define PI 3.14159265358979323846
#define MIN(a,b) ((a)<(b)?(a):(b))

int gsl_float_rotate(gsl_vector_float *hn, gsl_vector_float *he, float angle) {
/* rotates CCW in place seismic components hn, he which must form a RHS system by 
   angle (in deg)*/
  float c=cos(angle*PI/180.),s=sin(angle*PI/180.);
/*   printf("Rotate by angle %f c=%f s=%f\n",angle,c,s); */
/*   printf("Pre-rotation 0th element: %g %g\n",gsl_vector_float_get(hn,0),gsl_vector_float_get(he,0)); */
  gsl_blas_srot(he,hn,c,s);
/*   printf("Post-rotation 0th element: %g %g\n",gsl_vector_float_get(hn,0),gsl_vector_float_get(he,0)); */
  return(0);
}

void gsl_rotcr(gsl_vector *cor_mat[2][2], double azimuth) {
  /*% gsl_rotcr rotatin place2 D-cross-correlation matrix
    % crot=rotcr(C,azimuth) r otates cross-correlation matrix 
    % C=[ n1n1 e1n2 ; n1e2 e2e2 ] to direction azimuth (clockwise from N)
    NB Definition of angle is the opposite of gsl_float_rotate */ 
  double c=cos(azimuth*PI/180.),s=sin(azimuth*PI/180.);

  /* Left side mulltiplication with rotation matrix */
  gsl_blas_drot(cor_mat[0][0],cor_mat[1][0],c,s);
  gsl_blas_drot(cor_mat[0][1],cor_mat[1][1],c,s);
  /* Right side mulltiplication with rotation matrix transposed */
  gsl_blas_drot(cor_mat[0][0],cor_mat[0][1],c,s);
  gsl_blas_drot(cor_mat[1][0],cor_mat[1][1],c,s);
}
  


double gsl_df_normsqr(gsl_vector_float *a) {
  /* calculates Squared Norm=a_i a_i (summation convention) of vector */
  int i;
  double sum=0,b;
  for(i=0;i<a->size;i++) {
    b=(double)a->data[i*a->stride];
    sum+=b*b;
  }
  return(sum);
}

double gsl_df_dotprod(gsl_vector_float *a, gsl_vector_float *b) {
  /* calculates dot product=a_i b_i (summation convention) of vectors a,b */
  double result;
  gsl_blas_dsdot (a,b,&result);
  return(result);
}

void gsl_df_corrshift(gsl_vector *c, gsl_vector_float *x1, long beg1, gsl_vector_float *x2, long beg2, long len, long maxlag, int mode) {
  /*  calculates the modified cross-correlation, i.e. the dot product of vectors x1(b1-i/2:b1+len-1) and x2(b2+i/2:b2+len-1)
 For indices i between -maxlag and maxlag. x(j) for j<1 or j> length(x) is assumed to
 be zero.
 The result will be identical to the cross-correlation(a,b,maxlag) if x1,x2 are 0 outside the
 interval [bi:bi+len-1].

 mode: 
       GSL_CORRSHIFT_PARTITION 0     partition shift between both vectors
       GSL_CORRSHIFT_ONE 1           only shift the first sequence 
       GSL_CORRSHIFT_TWO 2           only shift the second sequence

 Output:
 c     the modified cross-correlation. C must be allocated outside the routine with
       a size of at least 2*maxlag+1  (if it is larger elements beyoun 2*maxlag remain
       unchanged.
 */
  gsl_vector_float_view vue1,vue2;
  long i, b1,b2,l;

  for(i=-maxlag;i<=maxlag;i++) {
    l=len;
    switch(mode) {
    case GSL_CORRSHIFT_PARTITION:
      b1=beg1+i/2-i;
      b2=beg2+i/2;
      break;
    case GSL_CORRSHIFT_ONE:
      b1=beg1-i;
      b2=beg2;
      break;
    case GSL_CORRSHIFT_TWO:
      b1=beg1;
      b2=beg2+i;
      break;
    }
    /* make sure bounds are within physical array.  If not truncate summations */
    if (b1<0) {
      b2 -= b1;  /* If part of first series is missing, then the corresponding parts
		    of the second would have been multiplied by two and also should not 
		    have been used */
      l -= b1;
      b1 =0;
    } 
    if (b2<0) {
      b1 -= b2;
      l  -= b2;
      b2 =0;
    }
    if(b1+l > x1->size) 
      l = x1->size - b1;
    if(b2+l > x2->size) 
      l = x2->size - b2;
    
/*     if (l!=len) printf("CORRSHIFT: adjusted parameters: i=%d b1 %d b2 %d l %d len %d beg1 %d len1 %d beg2 %d len2 %d\n",i,b1,b2,l,len,beg1,x1->size,beg2,x2->size); */
    
    vue1=gsl_vector_float_subvector(x1,b1,l);
    vue2=gsl_vector_float_subvector(x2,b2,l);
    gsl_blas_dsdot (&vue1.vector,&vue2.vector,gsl_vector_ptr(c,maxlag+i));  /* dot product of two sub-vectors i s the cross-correlation at that lag U*/
  }
}

void gsl_float_write_timeseries(FILE *fid,gsl_vector_float *data, float beg, float delta) {
  /* writes out a gsl_vector_float as xy file with  x being determined from time of first sample, beg,
     and sampling interval, delta */
  long i;
  for (i=0; i<data->size;i++)
    fprintf(fid,"%-20g %-20g\n",beg+i*delta,gsl_vector_float_get(data,i));
}

void gsl_float_write_pmp(FILE *fid,gsl_vector_float *data_x, gsl_vector_float *data_y) {
  /* writes out a two gsl_vector_float as parametric plot (e.g. for particle motion diagram) */
  long i,imax=MIN(data_y->size,data_x->size);
  for (i=0; i<imax;i++)
    fprintf(fid,"%-20g %-20g\n",gsl_vector_float_get(data_x,i), gsl_vector_float_get(data_y,i));
}

double gsl_df_polarisation(gsl_vector_float *data_n, gsl_vector_float *data_e, double *lin) {
  /* determine polarisation of a seismic signal */
  /* Return value: inferred polarisation in clockwise degree from N (0<pol<180 deg) */
  /* *lin:  Linearity   = 1-e2/e1 (e1>e2 eigenvalues of correlation matrix)
     lin = 1 for perfectly linear particle motion
     lin = 0 for circular motion */
  gsl_matrix *m2_cov=gsl_matrix_alloc(2,2);
  gsl_eigen_symmv_workspace *ws_eigen=gsl_eigen_symmv_alloc(2);
  gsl_matrix *m2_evec=gsl_matrix_alloc(2,2);
  gsl_vector *eval=gsl_vector_alloc(2);
  double pol;

  gsl_matrix_set(m2_cov,0,0,gsl_df_dotprod(data_n,data_n));
  gsl_matrix_set(m2_cov,0,1,gsl_df_dotprod(data_n,data_e));
  gsl_matrix_set(m2_cov,1,0,gsl_matrix_get(m2_cov,0,1));
  gsl_matrix_set(m2_cov,1,1,gsl_df_dotprod(data_e,data_e));

  gsl_eigen_symmv(m2_cov, eval, m2_evec, ws_eigen);              /* Calculate eigenvalues and vectors */
  gsl_eigen_symmv_sort(eval, m2_evec, GSL_EIGEN_SORT_VAL_ASC);   /* Sort them */
  
  *lin=1-gsl_vector_get(eval,0)/gsl_vector_get(eval,1);
  pol=atan2(gsl_matrix_get(m2_evec,1,1),gsl_matrix_get(m2_evec,0,1));
  return(pol*180/PI);
}
  

long nxtpwr2( long n) {
   long r; 
  for(r=1 ; r<n ; r <<= 1 ) { }
  return(r);
}

