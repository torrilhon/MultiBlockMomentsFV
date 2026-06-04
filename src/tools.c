#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#include "database_range.h"
#include "defs.h"
#include "tools.h"

// some tools for convenience, command line options and time measurements

struct timeval time0;


int exists_argument(int argc,const char **argv,const char *keystr)
{
  int i;

  for (i=1;i<argc;i++) if (strstr(argv[i],keystr)) return 1;
  return 0;
}

float get_float_argument(int argc,const char **argv,const char *keystr)
{
  int i;
  float result=0.0;

  for (i=1;i<argc;i++) {
    if (strstr(argv[i],keystr)) {
      result=atof(argv[i+1]);
      break;
    }
  }
  if (i == argc) {
    printf("Error: get_float_argument(%s)\n",keystr);
    return (float) -999999999;
  }
  else return result;
}

int get_string_argument(int argc,const char **argv,const char *keystr,char *result)
{
  int i;

  for (i=1;i<argc;i++) {
    if (strstr(argv[i],keystr)) {
      sprintf(result,"%s",argv[i+1]);
      break;
    }
  }
  if (i == argc) {
    fprintf(stderr,"Error: get_string_argument(%s)\n",keystr);
    return -1;
  }
  else return 0;
}




inline unsigned int GetmSec()
{
  struct timeval time;

  gettimeofday(&time, NULL );
  
  return( (unsigned int)(1.0e3*(time.tv_sec-time0.tv_sec)+1.0e-3*(time.tv_usec-time0.tv_usec)) );
}
  
void TickTimer()
{
  struct timeval time;

  gettimeofday(&time, NULL );
  
  if( time.tv_usec > time0.tv_usec ) {
    time0.tv_sec = time.tv_sec-time0.tv_sec;
    time0.tv_usec = time.tv_usec-time0.tv_usec;
  } else {
    time0.tv_sec = time.tv_sec-time0.tv_sec-1;
    time0.tv_usec = time0.tv_usec-time.tv_usec;
  };
};


int FullGridSize( DB_Block *B, int bMax )
{
  int nb, sum = 0;

  for( nb=0; nb<bMax; nb++ )
    sum += (B[nb].range.rbt.j-B[nb].range.lfb.j-2*MBC+1)*(B[nb].range.rbt.i-B[nb].range.lfb.i-2*MBC+1); 
  
  return( sum );
};

Vector GetPoint( double i, double j, double k, DB_Block B )
{
  Vector P;
  P.x = B.x+(i-B.range.lfb.i)*B.dx;
  P.y = B.y+(j-B.range.lfb.j)*B.dy;
  P.z = B.z;
  return( P );
};

