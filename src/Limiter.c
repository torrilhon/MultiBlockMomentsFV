//*********************************************************************************************
//                     ( slope limiter for second-order reconstruction )
//
//  The limiter variant is chosen at runtime via the global 'Limiter', which is set
//  from the block-file key 'limiter' or the -limiter flag (see LimiterCode below).
//
//  (c) Manuel Torrilhon, ETH Zurich, 2010 - RWTH Aachen 2026
//*********************************************************************************************

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "database_range.h"
#include "defs.h"
#include "tools.h"
#include "Limiter.h"


// limiter selection. Default: harmonic mean.
int Limiter = LIM_HARMONIC;

static double eps = 1.0e-10;


//********************** limiter functions for derivative reconstruction **********************
inline static double minmod( double d1, double d2 )
{
  return( 0.5*(sgn(d1)+sgn(d2))*min(fabs(d1),fabs(d2)) );
}

inline static double superbee( double theta )
{
  if( theta < 0.0 ) return( 0.0 );
  if( theta < 0.5 ) return( 2*theta );
  if( theta < 1.0 ) return( 1.0 );
  if( theta < 2.0 ) return( theta );
  return( 2.0 );
}

// slope limiter applied to the pair of one-sided differences (d1,d2).
// The variant is chosen at runtime via the global 'Limiter'.
double Limit1( double d1, double d2 )
{
  double w1, w2;
  switch( Limiter ) {
    case LIM_NONE:      // first order (no reconstruction)
      return( 0.0 );
    case LIM_AVERAGE:   // unlimited central average (second order, oscillatory)
      return( 0.5*(d1+d2) );
    case LIM_MINMOD:    // minmod (most diffusive TVD limiter)
      return( minmod(d1,d2) );
    case LIM_SUPERBEE:  // superbee (least diffusive TVD limiter)
      return( d1*superbee(d2/d1) );
    case LIM_HARMONIC:  // harmonic mean (van-Leer-type, smooth)
    default:
      w1 = fabs(d1)+eps; w2 = fabs(d2)+eps;
      return( (w1*d2+w2*d1)/(w1+w2) );
  };
}


//*********************************************************************************************
// map a limiter name (as used in the block file / -limiter flag) to its code,
// and back. Unknown names fall back to the harmonic-mean default.
//*********************************************************************************************
int LimiterCode( const char *name )
{
  if( strcmp(name,"none")     == 0 ) return( LIM_NONE );
  if( strcmp(name,"average")  == 0 ) return( LIM_AVERAGE );
  if( strcmp(name,"minmod")   == 0 ) return( LIM_MINMOD );
  if( strcmp(name,"superbee") == 0 ) return( LIM_SUPERBEE );
  if( strcmp(name,"harmonic") == 0 ) return( LIM_HARMONIC );
  fprintf( stderr, "   WARNING: unknown limiter '%s', using 'harmonic'.\n", name );
  return( LIM_HARMONIC );
}

const char *LimiterName( int code )
{
  switch( code ) {
    case LIM_NONE:     return( "none" );
    case LIM_AVERAGE:  return( "average" );
    case LIM_MINMOD:   return( "minmod" );
    case LIM_SUPERBEE: return( "superbee" );
    default:           return( "harmonic" );
  };
}
