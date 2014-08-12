#include <gsl/gsl_vector_float.h>
#include <math.h>
#include "sac_help.h"

#define F_EQ(a,b) (2.*fabs((a)-(b))/((a)+(b)+10*TOLERANCE)<TOLERANCE )

#define TOLERANCE 1e-5
#define VRBD(command) { if(1) { command  ; fflush(stdout); }}

static void sac_help_abort_msg(char *msg) {
  fprintf(stderr,"%s\n ABORT \n",msg);
  exit(10);
}


int check_consistency(sachdr *hdr1, sachdr *hdr2, unsigned int mode) {
/* checks whether two sac files are consistent. 

  Input:
  hdr1,hdr2  pointer to header structure of the two files to compare

  Flags:
  mode       What should I check.  Several modes can be combined by bit-OR
      CONSISTENCY_VERBOSE    report any inconsistency found to stderr;
      CONSISTENCY_EVENT      check reference time is the same, and evla,evlo, evdp are the same
      CONSISTENCY_STATION    check kstnm and knetwk and khole are  the same
      CONISTENCY_BEGIN       check that begin is the same
      Always check:          delta identical

  Returns: 0: SAC files are consistent
           1: SAC files are inconsistent
*/
  int verb=mode & CONSISTENCY_VERBOSE;
  if (!F_EQ(hdr1->delta,hdr2->delta)) {
/*    if (hdr1->delta != hdr2->delta ) { */
    if (verb) 
      fprintf(stderr,"Inconsistency in sampling rates:  %f vs %f\n",hdr1->delta,hdr2->delta);
    return(1);
  }
  if (mode & CONSISTENCY_EVENT) {
    if ( !F_EQ(hdr1->evla,hdr2->evla)
	 || !F_EQ(hdr1->evlo,hdr2->evlo)
	 || !F_EQ(hdr1->evdp,hdr2->evdp) ) {
      if (verb) 
	fprintf(stderr,"Inconsistency in event location:  (%f,%f,%f) vs (%f,%f,%f)\n",
		hdr1->evla,hdr1->evlo,hdr1->evdp,
		hdr2->evla,hdr2->evlo,hdr2->evdp);
      return(1);
    }
    if ( hdr1->nzyear != hdr2->nzyear 
	 || hdr1->nzjday != hdr2->nzjday 
	 || hdr1->nzhour != hdr2->nzhour
	 || hdr1->nzmin  != hdr2->nzmin 
	 || hdr1->nzsec  != hdr2->nzsec ) {
      if (verb)
	fprintf(stderr,"Inconsistency in reference time\n");
      return(1);
    }
  }
  if (mode & CONSISTENCY_STATION) {
    if ( strncmp(hdr1->knetwk,hdr2->knetwk,8) 
	 || strncmp(hdr1->kstnm,hdr2->kstnm,8) 
	 || strncmp(hdr1->khole,hdr2->khole,8) ) {
      if (verb)
	fprintf(stderr,"Inconsistency in station information (%.8s,%.8s,%.8s) vs (%.8s,%.8s,%.8s)\n",
		hdr1->knetwk,hdr1->kstnm,hdr1->khole,
		hdr2->knetwk,hdr2->kstnm,hdr2->khole);
      return(1);
    }
  }
  if (mode & CONSISTENCY_BEGIN) {
    if ( !F_EQ(hdr1->b,hdr2->b)) {
      if (verb) 
	fprintf(stderr,"Inconsistency in seismogram begin time:  %f vs %f\n",
		hdr1->b,hdr2->b);
      return(1);
    }
  }

  return(0);
}

void make_rhs(sachdr **hdr1, gsl_vector_float **data1, sachdr **hdr2, gsl_vector_float **data2) {
  /* checks that horizontal azimuths differ by 90 degree and exit with an error message if they are not
  .  If the second component is not to the right of the first component, exchange the components */
  sachdr *dummy;
  gsl_vector_float *vector_dummy;
  if (F_EQ(fmod((*hdr1)->cmpaz-(*hdr2)->cmpaz+180,360)-180,90)) {
    /* components form LHS system, need to mirror one of them */
    warn_msg("Horizontal components form LHS system. Exchanging components.\n");
    dummy=*hdr1;
    *hdr1=(*hdr2);
    (*hdr2)=dummy;
    vector_dummy=*data1;
    *data1=*data2;
    *data2=vector_dummy;
  }
  else if (!F_EQ(fmod((*hdr2)->cmpaz-(*hdr1)->cmpaz+180,360)-180,90)) {
    sac_help_abort_msg("Horizontal components not perpendicular");
  }
}

void read_seis_file(char *fname, sachdr **hdr, gsl_vector_float **data) {
/* read header and data of SAC file. Function allocates spaces for header and data
Input: char *fname       File name of sac file

   Ouput: sachdr *hdr      Sac header structure
          gsl_vector_float *data Vector containing seismic data


   Note: this function reads the sac header directly without using
         SAC library functions.  It might not work with future versions
         of SAC files 
*/
  FILE *in;
  float dummy;
  if((in=fopen(fname,"rb"))==NULL) {
    fprintf(stderr,"Error opening SAC file %s.",fname); sac_help_abort_msg("");
  }

  *hdr=(sachdr *)malloc(sizeof(sachdr));

  if (fread(*hdr,sizeof(sachdr),1,in)!=1) {
    fprintf(stderr,"Error reading SAC file %s.",fname); sac_help_abort_msg("");
  }
  /* Check that this is really a SAC file */
  if((*hdr)->internal2 != -12345. || (*hdr)->internal3 != -12345. || (*hdr)->nvhdr > 6) {
    fprintf(stderr,"%s does not appear to be a SAC file or has wrong byte order.",fname);
    sac_help_abort_msg("");
  }
  /*  VRBD(printf("DEBUG: Reading sacfile %s length %d, year %d,nvhdr %d \n",fname,(*hdr)->npts,(*hdr)->nzyear,fname,(*hdr)->nvhdr)); */
  *data=gsl_vector_float_alloc((*hdr)->npts);  /* NOTE: gsl error routine takes care off files which are too short */
  gsl_vector_float_fread(in,*data);
  fread(&dummy,sizeof(float),1,in);
  if (!feof(in)) {
    fprintf(stderr,"File length too long for NPTS value %d in sac file %s",(*hdr)->npts,fname);
    gsl_vector_float_free(*data);
    sac_help_abort_msg("");
  }
}

void write_seis_file(char *fname, sachdr *hdr, gsl_vector_float *data) {
/* read header and data of SAC file.  
Input: char *fname       File name of sac file

       sachdr *hdr      Sac header structure
       gsl_vector_float *data Vector containing seismic data


   Note: this function writes the sac header directly without using
         SAC library functions.  It might not work with future versions
         of SAC files 
*/
  FILE *out;
  float dummy;
  if((out=fopen(fname,"wb"))==NULL) {
    fprintf(stderr,"Error opening SAC file %s.",fname); sac_help_abort_msg("");
  }

  if (fwrite(hdr,sizeof(sachdr),1,out)!=1) {
    fprintf(stderr,"Error writing SAC file %s.",fname); sac_help_abort_msg("");
  }
  if (gsl_vector_float_fwrite(out,data) != 0) {
    fprintf(stderr,"Error writing SAC file %s.",fname); sac_help_abort_msg("");
  }
}

