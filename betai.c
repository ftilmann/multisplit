#include <gsl/gsl_sf_gamma.h>
#include <stdio.h>
#include <stdlib.h>

/* wrapper routine for gsl_betai function using NR syntac */
double betai(double a, double b, double x)
{
  if (x < 0.0 || x > 1.0) {
    fprintf(stderr,"Numerical error: Bad x (%f) in routine betai\n",x);
    exit(10);
  }
  return gsl_sf_beta_inc(a,b,x);
}
