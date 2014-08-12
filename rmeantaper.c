#include <stdio.h>
#include <math.h>

void rmean_and_taper(float data[], long int len, float del, int rmean, float taperlen) {
  long i;
  float mean,fac;
  long talen;

  /*  fprintf(stderr,"rmean_and_taper %d %f\n",rmean,taperlen); */
  if (rmean) {
    mean=0;
    for (i=0; i<len; i++)
      mean+=data[i]/len;
/*     fprintf(stderr,"Remove Mean: %f Length: %d\n",mean,len); */
    for (i=0; i<len; i++)
      data[i]-=mean;
  }
  if (taperlen>0.0) {
/*     fprintf(stderr,"Taper\n"); */
    talen=(unsigned long int)(taperlen/del);
    if (talen>len/2) 
      talen=len/2;
    for (i=0; i<talen; i++) {
      fac=sin((i+1)*3.1415926535897/(2.*(float)(talen+1)));
      data[i]*=fac;
      data[len-i-1]*=fac;
    }
  }
}    
