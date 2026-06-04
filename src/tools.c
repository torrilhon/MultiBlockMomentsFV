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


// Find flag `keystr` and return the index of its VALUE token (i+1), or -1
// if the flag is absent or is the last token with no value following it.
// Exact match on the flag (so "-out" does not match "-output") and the
// bounds check on i+1 avoid the out-of-range argv read.
static int argument_value_index(int argc,const char **argv,const char *keystr)
{
  int i;

  for (i=1;i<argc;i++)
    if (strcmp(argv[i],keystr) == 0)
      return (i+1 < argc) ? i+1 : -1;
  return -1;
}

int exists_argument(int argc,const char **argv,const char *keystr)
{
  int i;

  for (i=1;i<argc;i++) if (strcmp(argv[i],keystr) == 0) return 1;
  return 0;
}

float get_float_argument(int argc,const char **argv,const char *keystr)
{
  int v = argument_value_index(argc,argv,keystr);
  if (v < 0) {
    printf("Error: get_float_argument(%s)\n",keystr);
    return (float) -999999999;
  }
  return atof(argv[v]);
}

int get_string_argument(int argc,const char **argv,const char *keystr,char *result)
{
  int v = argument_value_index(argc,argv,keystr);
  if (v < 0) {
    fprintf(stderr,"Error: get_string_argument(%s)\n",keystr);
    return -1;
  }
  sprintf(result,"%s",argv[v]);
  return 0;
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

