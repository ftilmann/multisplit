/* split_cor */
/* Apply or undo splitting specified by fast direction and splitting
   time to two horizontal component SAC files */
/* Note length of SAC file and B header will change as both fast and 
   slow samples must be available for the output file */
/* Author: F Tilmann */
/* Contact: tilmann|a|gfz-potsdam.de */
/* (C) 2004 F Tilmann */
/* This source code is released under the GNU public license */

/* Code uses sac, gsl and gslblas library and needs GMT executables to be installed */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <gsl/gsl_vector_float.h>

#include "gsl_seis.h"
#include "sac_help.h" 

#define MAXLEN 4194304

typedef struct {
  /* Input data */
  sachdr *hdr_hor1;
  gsl_vector_float *data_hor1; 
  char *file1;
  sachdr *hdr_hor2;
  gsl_vector_float *data_hor2;
  char *file2;
  float fast;
  float time;
  char pstr[32];
  char astr[32];
} params;

int verbose=0;
#define VRB(command) { if(verbose) { command  ; fflush(stdout); }}


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
                         should be measured.  The header variables CMPAZ must be set.  Note that\n\
                         the whole time series is splitting-corrected even though this might\n\
                         not always be appropriate\n\
\n\
-split fast time         Give fast direction in degree clockwise from N, and splitting time in s\n\
\n\
Choose at least one of:\n\
\n\
-pre pstr                Generate output file names by prepending pstr to filenames\n\
\n\
-ap astr                Generate output file names by appending pstr to filenames\n\
\n\
Optional modifiers:\n\
\n\
-forward                 Apply splitting to waveform\n\
                         (Default is to remove specified splitting from waveform)\n\
\n\
Status options:\n\
\n\
-v                       Verbose output (for debugging)\n\
-h                       Print out this help text\n\
");                         
  exit(10);
}

void parse(int argc, char **argv, params *par) {
  int iarg;
  char *dummy;
  int read_data=0, invert=1;
  /* Required arguments: initialise and check later */
  par->fast=-999;
  par->time=-999;  
  par->pstr[0]='\0';
  par->astr[0]='\0';
  par->pstr[31]='\0';
  par->astr[31]='\0';

  iarg=0;
  if (argc==1) 
    usage("split_cor");
  while(++iarg<argc) {
    if(argv[iarg][0]!='-') {
      fprintf(stderr,"%s ",argv[iarg]);
      abort_msg("is not an option. This error can also be due to the wrong number of arguments for a previous option"); 
    }
    if(!strcasecmp(argv[iarg],"-data")) {
      read_data=TRUE;
      if ( iarg+2>=argc ) 
	abort_msg("-data option must be followed by two sac-filenames (horizontal components).");
      read_seis_file(argv[++iarg],&(par->hdr_hor1),&(par->data_hor1));
      par->file1=strdup(argv[iarg]);
      read_seis_file(argv[++iarg],&(par->hdr_hor2),&(par->data_hor2));
      par->file2=strdup(argv[iarg]);
    }
    else if(!strncasecmp(argv[iarg],"-pre",4)) {
      iarg++;
      strncpy(par->pstr,argv[iarg],31);
    }
    else if(!strncasecmp(argv[iarg],"-app",4)) {
      iarg++;
      strncpy(par->astr,argv[iarg],31);
    }
    else if(!strncasecmp(argv[iarg],"-spl",3)) {
      par->fast=atof(argv[++iarg]);
      par->time=atof(argv[++iarg]);
    }
    else if(!strncasecmp(argv[iarg],"-for",4)) {
      invert=0;
    }
    else if(!strcasecmp(argv[iarg],"-v") ) {
      verbose=1;
    }
    else if(!strcasecmp(argv[iarg],"-h") ) {
      usage("split_cor");
    }
    else {
      fprintf(stderr,"%s ",argv[iarg]);
      abort_msg("is not a known option");
    }
  }
  if (!read_data)
    abort_msg("Need to specify -data option");
  if (par->fast==-999) 
    abort_msg("Need to specify fast direction and splitting time");
  if (strlen(par->pstr)==0 && strlen(par->astr)==0) 
    abort_msg("Need to define at least one of -pstr, -astr");
  if(invert==1) 
    par->time = -par->time;
}



int main(int argc, char **argv)
{
  params *par=(params *) malloc(sizeof(params));
  char output[256];
  long ishift;
  gsl_vector_float_view vue_fast,vue_slow;

  gsl_vector_float *hor1,*hor2;

  parse(argc,argv, par);

  /* Check that sac files are consistent */
  if(check_consistency(par->hdr_hor1,par->hdr_hor2, 
                       CONSISTENCY_STATION | CONSISTENCY_EVENT | CONSISTENCY_VERBOSE | CONSISTENCY_BEGIN )) 
    abort_msg("Data files inconsistent");     
  make_rhs(&par->hdr_hor1,&par->data_hor1,&par->hdr_hor2,&par->data_hor2);

  if (par->time < 0) {
    par->time = -par->time;
    par->fast = par->fast+90;
  }
  /* rotate to fast/slow */
  gsl_float_rotate(par->data_hor1,par->data_hor2, par->hdr_hor1->cmpaz - par->fast);
  ishift=(long)rint(par->time/par->hdr_hor1->delta);
  VRB(printf("Shift %f ishift %d New length %d\n",par->time,ishift,par->data_hor1->size-ishift));
  vue_fast=gsl_vector_float_subvector(par->data_hor1,ishift,par->data_hor1->size-ishift);
  vue_slow=gsl_vector_float_subvector(par->data_hor2,0,par->data_hor2->size-ishift);
  par->hdr_hor1->npts -= ishift;
  par->hdr_hor2->npts -= ishift;
  par->hdr_hor1->b += par->hdr_hor1->delta*ishift/2;
  par->hdr_hor2->b += par->hdr_hor2->delta*ishift/2;
  
  /* rotate back to CMPAZ (also mangles original data vectors */
  gsl_float_rotate(&vue_fast.vector,&vue_slow.vector, -par->hdr_hor1->cmpaz + par->fast);

  /* and write out */
  strcpy(output,par->pstr);
  strcat(output,par->file1);
  strcat(output,par->astr);
  write_seis_file(output,par->hdr_hor1,&vue_fast.vector);
  
  strcpy(output,par->pstr);
  strcat(output,par->file2);
  strcat(output,par->astr);
  write_seis_file(output,par->hdr_hor2,&vue_slow.vector);

}




