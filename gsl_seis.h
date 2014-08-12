#include <math.h>
#include <gsl/gsl_vector_float.h>
#include <gsl/gsl_vector_double.h>


int gsl_float_rotate(gsl_vector_float *hn, gsl_vector_float *he, float angle);
void gsl_float_write_pmp(FILE *fid,gsl_vector_float *data_x, gsl_vector_float *data_y);
void gsl_float_write_timeseries(FILE *fid,gsl_vector_float *data, float beg, float delta);

double gsl_df_polarisation(gsl_vector_float *data_n, gsl_vector_float *data_e, double *lin);
double gsl_df_normsqr(gsl_vector_float *a);
double gsl_df_dotprod(gsl_vector_float *a, gsl_vector_float *b);
void gsl_df_corrshift(gsl_vector *c, gsl_vector_float *x1, long beg1, gsl_vector_float *x2, long beg2, long len, long maxlag, int mode);

void gsl_rotcr(gsl_vector *cor_mat[2][2], double azimuth);

long nxtpwr2( long n) ;

/* modes gos gsl_df_corrshift */
#define GSL_CORRSHIFT_PARTITION 0  /* partition shift between both vectors */
#define GSL_CORRSHIFT_ONE 1           /* only shift the first sequence */
#define GSL_CORRSHIFT_TWO 2           /* only shift the second sequence */
