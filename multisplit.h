#include <gsl/gsl_vector_float.h>
#include "gsl_seis.h"
#include "sac.h"

#define SQR(a) ((a)*(a))
#define CUB(a) ((a)*(a)*(a))
#define ABS2(a,b) ((a)*(a)+(b)*(b))
#define ABS(real,im) (sqrt((real)*(real)+(im)*(im)))
#define ARG(real,im) (atan2((im),(real)))

#define F_EQ(a,b) (2.*fabs((a)-(b))/((a)+(b)+10*TOLERANCE)<TOLERANCE )

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

/* Note: my system at home does not appear to have the round function */
#define ROUND(a) (floor((a)+0.5))

/* PER(per,x) return a value v with range 0..per-1 such that v+n*per=x (n an integer). It only works for integers */
#define PER(per,x) ( x<0 ? (x + per*(1+(-x)/per)) : (x%per) )

#define PI 3.14159265358979323846
#define TOLERANCE 1e-5
#ifndef NAN
#define NAN (strtod("NAN",NULL))
#endif

#define CMAT_MEMCPY(dest,src) { \
  gsl_vector_memcpy(dest[0][0],src[0][0]); \
  gsl_vector_memcpy(dest[0][1],src[0][1]); \
  gsl_vector_memcpy(dest[1][0],src[1][0]); \
  gsl_vector_memcpy(dest[1][1],src[1][1]); }

#define CMAT_ALLOCATE(mat,length) { \
  mat[0][0]=gsl_vector_alloc(length); \
  mat[0][1]=gsl_vector_alloc(length); \
  mat[1][0]=gsl_vector_alloc(length); \
  mat[1][1]=gsl_vector_alloc(length);}

#define CMAT_FREE(mat) { \
  gsl_vector_free(mat[0][0]); \
  gsl_vector_free(mat[0][1]); \
  gsl_vector_free(mat[1][0]); \
  gsl_vector_free(mat[1][1]); }

#define TRUE 1
#define FALSE 0

typedef int logical;


/* Method specific parameters */
typedef struct {
  float maxshift;
  sachdr *hdr_ref1;
  gsl_vector_float *data_ref1; 

  sachdr *hdr_ref2;
  gsl_vector_float *data_ref2;
} correl_params;

/* Model specific parameters */
typedef struct {
  float faststep;
  float fastmin;
  float fastmax;
  float timemin;
  float timemax;
  float timestep;
} hor_split;

typedef struct {
  hor_split top; /* splitting parameters of top  layer */
  hor_split bot; /* splitting parameters of bottom layer */
} split_params;

typedef struct {
  char phase_start[9];
  char phase_end[9];
  float offset_start;
  float offset_end;
  float taper ;    /* if taper<0 then do not taper but move windows through data */
} analysis_window;

/* Global parameters */
typedef struct {
  /* Input data */
  sachdr *hdr_hor1;
  gsl_vector_float *data_hor1; 
  sachdr *hdr_hor2;
  gsl_vector_float *data_hor2;
  /* Method */
  int method ;
  union  {
    correl_params cor_par;
  } method_q;
  /* Model space definition */
  int model;
  union {
    split_params split_par;
  }  model_q;
  /* Window definition */
  analysis_window window;
  /* Modifiers */
  int make_grd ;    /* MAKE_GRD : create GMT grid files MAKE_GMT: create and execute full GMT scripts */
  float dof_s ;         /* degrees of freedom per second */
  char root[128] ;      /* Root for filename */  
} ms_params  ;

/* METHOD IDENTIFIERS */
#define MINEVALUE 1
#define MINTRANSVERSE 2
#define CONV 3
#define CORREL 4

/* MODEL IDENTIFIERS */
#define SINGLE_HOR_SPLIT 1
#define DOUBLE_HOR_SPLIT 2

/* make_grd identifieres */
#define MAKE_GRD 1
#define MAKE_GMT 2
#define MAKE_GMT5 2|4   /* note that MAKE_GMT5 implies MAKE_GMT (second bit set) */

/* Function protoptypes */
void abort_msg(char *);
void err_single_split_sks(ms_params *par, gsl_matrix *m_res_energy, gsl_vector_float *north, gsl_vector_float *east,long beg, long len, long mod_par, char *phase,  gsl_matrix *m_aux1, gsl_matrix *m_aux2, gsl_vector_float *ref_north, gsl_vector_float *ref_east, long refbeg);
float find_phase(sachdr *hdr,char *phase_name, char *phase_out);
gsl_vector_float *find_window(sachdr *hdr, gsl_vector_float *data, analysis_window *win, float maxlag, long *beg, long *len, char *phase);
double invfisher(double nu1,double nu2,double conf);
void make_event_name(char *string, sachdr *hdr, int mode);
FILE *open_for_write(char *root,char *extension);
void parse(int argc, char **argv, ms_params *par);
void single_split_correl(int method, hor_split *hsplit, float maxshift, gsl_vector_float *nspl, gsl_vector_float *espl, long beg, gsl_vector_float *ref_north, gsl_vector_float *ref_east, long ref_beg, long len, float delta, gsl_matrix *m_res_energy, gsl_matrix *m_delay, gsl_matrix *m_alpha);
void single_split_sks(int method, hor_split *hsplit, gsl_vector_float *north, gsl_vector_float *east, long beg, long len, float delta, float baz, gsl_matrix *res_energy, gsl_matrix *pol);
void usage(char *cmd);
void warn_msg(char *);


/* Modes for make_event_name */
#define EVN_YYJJJHHMM 1

/* External function prototypes */
void rmean_and_taper(float data[], long int len, float del, int rmean, float taperlen);
double gammln(double xx);
double betai(double a, double b, double x);
double betacf(double a, double b, double x);
