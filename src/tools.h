#include <sys/time.h>

static inline double max( double a, double b ) 
{	
  return( a > b ? a : b );
}

static inline double min( double a, double b ) 
{	
  return( a < b ? a : b );
}

static inline double sgn( double a ) 
{ 
  if( a == 0 ) return( 0 );
  return( a < 0 ? -1 : 1 );
}

int exists_argument(int argc,const char **argv,const char *keystr);
float get_float_argument(int argc,const char **argv,const char *keystr);
int get_string_argument(int argc,const char **argv,const char *keystr,char *result);


extern struct timeval time0;
unsigned int GetmSec();
void TickTimer();

int FullGridSize( DB_Block *B, int bMax );
Vector Get_Point( double i, double j, double k, DB_Block B );


