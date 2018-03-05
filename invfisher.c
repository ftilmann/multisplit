#include <stdio.h>
#include <stdlib.h>

#define MAX(a,b) ((a)>(b)?(a):(b))

double betai(double a, double b, double x);

double invfisher(double nu1,double nu2,double conf) {
  /* note this routine is transliterated from MATLAB function */
  /*  function e95=invfisher(nu1,nu2,conf)
      % invfisher calculates inverse Fisher distribution
      % e95=invfisher(nu1,nu2,x)
      % Perform bisection to find the the x confidence interval for
      % degrees of freedom nu1,nu1
      % x defaults to 0.95
      %
      % Taken from M.Bostock, T.Hearn's code erranal.m for splitting analysis
      % Generalised by F Tilmann
  */
  int jmax=40,j;
  double x1,x2,f,fmid,rtbis,dx,xmid;
  double xacc=0.01*nu2/(nu2+nu1*200.0);
  /*%
    % Bracket the root between values of F = 1 and 100. Note
    % that we are determining the f(1-a) confidence region
    % a=0.05. */
  x1=nu2/(nu2+nu1*10000.0); 
  f=(1.0-betai(nu2/2.0,nu1/2.0,x1))-conf; 
  x2=nu2/(nu2+nu1*1.0); 
  fmid=(1.0-betai(nu2/2.0,nu1/2.0,x2))-conf; 
  if (f*fmid >= 0) {
    fprintf(stderr,"ERROR invfisher: root must be bracketed");
    exit(10);
  }
  if (f < 0) {
    rtbis=x1;
    dx=x2-x1;
  } else {
    rtbis=x2;
    dx=x1-x2;
  }
  j=0;
  while (j < jmax) {
    dx=dx*0.5;
    xmid=rtbis+dx;
    fmid=(1.0-betai(nu2/2.0,nu1/2.0,xmid))-conf;
    if (fmid <= 0) 
      rtbis=xmid;
    j=j+1;
    if ( MAX(dx,-dx) < xacc || fmid == 0) 
      break; 
  }
  if (j>jmax) 
    fprintf(stderr,"WARNING invfisher: root finding did not converge. Error estimates likely to be meaningless");
  return (nu2/xmid-nu2)/nu1;
}
