#include <gsl/gsl_vector_float.h>
#include "sac.h"
typedef struct sac sachdr;

int make_abort(int hello);
int check_consistency(sachdr *hdr1, sachdr *hdr2, unsigned int mode);
void make_rhs(sachdr **hdr1, gsl_vector_float **data1, sachdr **hdr2, gsl_vector_float **data2);
void read_seis_file(char *fname, sachdr **hdr, gsl_vector_float **data);
void write_seis_file(char *fname, sachdr *hdr, gsl_vector_float *data);

/* Modes for check_consistency */
#define CONSISTENCY_VERBOSE  1<<1
#define CONSISTENCY_STATION  1<<2
#define CONSISTENCY_EVENT    1<<3
#define CONSISTENCY_BEGIN    1<<4
