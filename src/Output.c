//*********************************************************************************************
//                     ( Tecplot output of the solution )
//
//  Writes Tecplot POINT-format files. 'Convert' formats a double into a compact
//  filename token (e.g. 0.05 -> "0p05").
//
//  (c) Manuel Torrilhon, ETH Zurich, 2010 - RWTH Aachen 2026
//*********************************************************************************************

#define _GNU_SOURCE     // asprintf
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "database_range.h"
#include "defs.h"
#include "tools.h"
#include "Flux.h"
#include "Scenario.h"   // Kn, R13, Type
#include "Output.h"


char *Convert( double a, char *s )
{
  int i = fabs(a+1e-7);
  int j = fabs(10*(a-i+1e-7));
  int k = fabs(100*(a-i-0.1*j+1e-7));
  int l = fabs(1000*(a-i-0.1*j-0.01*k+1e-7));

  if( l != 0 ) snprintf( s, 20, "%dp%d%d%d", i, j, k, l );
  else if( k != 0 ) snprintf( s, 20, "%dp%d%d", i, j, k );
  else snprintf( s, 20, "%dp%d", i, j );
  return( s );
};

//*********************************************************************************************
// writes tecplot-files with multiple time results
// Add > 0 : appends data to existing tAll-file with ID=Add+1
// Add = 0 : Creates new file with time tXXpXX
// Add = -1 : Creates new file with tAll-file and ID=1
//*********************************************************************************************
void Output( DB_Block *B, int bMax, int n, double t, int Add )
{
  extern char OutDir[200];   // output directory set via -out ("" = current dir)
  int i, j, k = 1, nb, dst_idx, nGhost = 1, lineOut;
  double Prim[N];
  char s1[20], s2[20], solv[1024], prefix[224];
  char *name = NULL;
  Vector P;
  FILE *fp;

  // prepend output directory (with trailing slash) if -out was given
  if( strlen(OutDir) > 0 ) snprintf( prefix, sizeof prefix, "%s/", OutDir );
  else prefix[0] = '\0';

  if( R13 < 1.0 )
    snprintf( solv, sizeof solv, "%s%s", "Grad13", Type );
  else
    snprintf( solv, sizeof solv, "%s%s", "Reg13", Type );

  if( !Add ) Convert(t,s2);
  else snprintf( s2, sizeof s2, "All" );

  if( asprintf(&name, "%s%s_nst%d_Kn%s_t%s.dat", prefix, solv,
    FullGridSize( B, bMax ), Convert(Kn,s1), s2 ) == -1 ) name = NULL;
  
  if( Add <= 0 ) {
    Add = 0;
    fp = fopen( name, "w" );
    if( fp == NULL ) { perror( name ); if(name) free(name); return; }
    fprintf( fp, "TITLE = \"13 Moment Equations in 2D\"\n" ); 
    fprintf( fp, "VARIABLES = \"x\",\"y\",\"rho\",\"vx\",\"vy\",\"T\",\"p\","); 
    fprintf( fp, "\"Pxx\",\"Pxy\",\"Pyy\",\"qx\",\"qy\"\n"); 
  } else {
    fp = fopen( name, "a+" );
    if( fp == NULL ) { perror( name ); if(name) free(name); return; }
  };
  
  if( Add == 0 ) printf( "\n  Create File: %s...", name );
  if( Add > 0 )  printf( "\n  Append (%d) File: %s...", Add, name );
  fflush(stdout);

  for( nb=0; nb<bMax; nb++ ) {
    lineOut = B[nb].range.rbt.i - B[nb].range.lfb.i + 1;
    fprintf( fp, "ZONE T=\"t = %5.3f, Block %d\", I=%d, J=%d, F=POINT", t, nb+1, 
          B[nb].range.rbt.j-B[nb].range.lfb.j-2*(MBC-nGhost)+1, 
          B[nb].range.rbt.i-B[nb].range.lfb.i-2*(MBC-nGhost)+1 ); 
    fprintf( fp, ", SOLUTIONTIME=%5.3f, STRANDID=%d\n", t, Add+1 );
    
    for( i=B[nb].range.lfb.i+MBC-nGhost; i<=B[nb].range.rbt.i-MBC+nGhost; i++ ) {
      for( j=B[nb].range.lfb.j+MBC-nGhost; j<=B[nb].range.rbt.j-MBC+nGhost; j++ ) {
        P = GetPoint(i,j,k,B[nb]);
        dst_idx = (j - B[nb].range.lfb.j)*lineOut + (i - B[nb].range.lfb.i);
        Primitive( B[nb].layer[L_DATA][dst_idx].U, Prim );

        fprintf( fp, "%20.8f", P.x );
        fprintf( fp, "%20.8f", P.y );
        fprintf( fp, "%20.8f", Prim[0] );
        fprintf( fp, "%20.8f", Prim[1] );
        fprintf( fp, "%20.8f", Prim[2] );
        fprintf( fp, "%20.8f", Prim[3]/Prim[0] );
        fprintf( fp, "%20.8f", Prim[3] );
        fprintf( fp, "%20.8f", Prim[4] );
        fprintf( fp, "%20.8f", Prim[5] );
        fprintf( fp, "%20.8f", Prim[6] );
        fprintf( fp, "%20.8f", Prim[7] );
        fprintf( fp, "%20.8f", Prim[8] );
        fprintf( fp, "\n" );
      };
    };
  };
  
  fclose( fp );
  if( name ) free(name);
  printf( " done\n\n" );

};
