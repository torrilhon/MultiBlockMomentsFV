//*********************************************************************************************
//                     ( binary output of the solution )
//
//  Writes a two-file sidecar result per run, sharing one basename:
//
//    <name>.hdr   pure ASCII run manifest -- human-readable / greppable, holds the
//                 resolved run parameters plus one ZONE record per (block,time)
//                 carrying I, J, SOLUTIONTIME, STRANDID and the byte OFFSET/NBYTES
//                 of that zone's payload inside the .dat file.
//    <name>.dat   pure binary -- nothing but contiguous Float32 (default) or
//                 Float64 payloads for all blocks and all appended times,
//                 concatenated back-to-back.
//
//  The 2 grid coordinates x,y are NOT stored: the grids are uniform Cartesian,
//  so the .hdr records each block's geometry (BX,BY,DX,DY,LFB_I,LFB_J,ISTART,
//  JSTART) and the reader reconstructs coordinates exactly from GetPoint's
//  linear map. Only the 10 physical fields are written per point.
//
//  Layout per zone (variable-fastest, Julia-friendly):
//      rho,vx,vy,T,p,Pxx,Pxy,Pyy,qx,qy  for point 1,
//      rho,vx,vy,...                     for point 2, ...
//  so a zone reads back as  reshape(read!(io, Vector{T}(undef, NVARS*I*J)),
//  NVARS, I*J)  with columns == points and T = Float32/Float64 per PRECISION.
//
//  Output precision defaults to Float32 (single); pass -double for Float64.
//
//  'Convert' formats a double into a compact filename token (e.g. 0.05 -> "0p05").
//
//  Add > 0 : appends a new time instance to the existing tAll files (STRANDID = Add+1)
//  Add = 0 : creates new files for a single time tXXpXX
//  Add = -1: creates new tAll files with STRANDID = 1
//
//  (c) Manuel Torrilhon, ETH Zurich, 2010 - RWTH Aachen 2026
//*********************************************************************************************

#define _GNU_SOURCE     // asprintf
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "database_range.h"
#include "defs.h"
#include "tools.h"
#include "Flux.h"
#include "Scenario.h"   // Kn, R13, Type
#include "Mesh.h"       // Refine
#include "Output.h"

#ifndef GIT_COMMIT
#define GIT_COMMIT "unknown"
#endif

#define NVARS_OUT 10               // 10 physical fields per point (x,y reconstructed)
#define FORMAT_VERSION 1

int OutFloat32 = 1;                 // 1 = single-precision output (default), 0 = double
                                   // ('extern' in Output.h, set in Main from -double)

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
// writes the ASCII manifest preamble (provenance + [RUN] + [VARIABLES] + [ZONES] tag).
// Called only when a fresh .hdr is created (Add <= 0).
//*********************************************************************************************
static void WriteManifest( FILE *fh, const char *datname, const char *solv,
                           int gridSize, double dt, double tEnd )
{
  char tbuf[64];
  time_t now = time(NULL);
  struct tm *utc = gmtime(&now);
  strftime( tbuf, sizeof tbuf, "%Y-%m-%dT%H:%M:%SZ", utc );

  fprintf( fh, "# MultiBlockMomentsFV binary output header\n" );
  fprintf( fh, "# FORMAT_VERSION = %d\n", FORMAT_VERSION );
  fprintf( fh, "# CREATED_UTC = %s\n", tbuf );
  fprintf( fh, "# GIT_COMMIT = %s\n\n", GIT_COMMIT );

  fprintf( fh, "[RUN]\n" );
  fprintf( fh, "solver = %s\n", solv );
  fprintf( fh, "Kn = %.10g\n", Kn );
  fprintf( fh, "R13 = %.10g\n", R13 );
  fprintf( fh, "Type = %s\n", Type );
  fprintf( fh, "FullGridSize = %d\n", gridSize );
  fprintf( fh, "refine = %d\n", Refine );
  fprintf( fh, "dt = %.10g\n", dt );
  fprintf( fh, "t = %.10g\n\n", tEnd );

  fprintf( fh, "[BINARY]\n" );
  fprintf( fh, "BINARY_FILE = %s\n", datname );
  fprintf( fh, "NVARS = %d\n", NVARS_OUT );
  fprintf( fh, "ORDER = variable-fastest\n" );
  fprintf( fh, "PRECISION = %s\n", OutFloat32 ? "Float32" : "Float64" );
  fprintf( fh, "ENDIAN = little\n\n" );

  // x,y are reconstructed by the reader from the per-zone geometry (uniform
  // Cartesian cell centres); only the 10 physical fields are stored in .dat.
  fprintf( fh, "[VARIABLES]\n" );
  fprintf( fh, "rho vx vy T p Pxx Pxy Pyy qx qy\n\n" );

  fprintf( fh, "[ZONES]\n" );
  fprintf( fh, "# block I J NVARS BX BY DX DY LFB_I LFB_J ISTART JSTART NGHOST SOLUTIONTIME STRANDID OFFSET NBYTES\n" );
};


//*********************************************************************************************
// writes binary output of the solution as a .hdr (ASCII manifest) + .dat (binary) pair.
//*********************************************************************************************
void Output( DB_Block *B, int bMax, int n, double t, double dt, double tEnd, int Add )
{
  extern char OutDir[200];   // output directory set via -out ("" = current dir)
  int i, j, nb, dst_idx, nGhost = 1, lineOut, ni, nj, i0, j0, p, v, m;
  size_t elemSize = OutFloat32 ? sizeof(float) : sizeof(double);
  long offset;
  long long nbytes;
  double Prim[N], vals[NVARS_OUT];
  char s1[20], s2[20], solv[1024], prefix[224];
  char *base = NULL, *hdrname = NULL, *datname = NULL;
  const char *datbase;            // .dat basename without directory, for the manifest
  unsigned char *buf = NULL;      // raw staging buffer (Float32 or Float64 per OutFloat32)
  FILE *fh, *fd;
  (void)n;

  // prepend output directory (with trailing slash) if -out was given
  if( strlen(OutDir) > 0 ) snprintf( prefix, sizeof prefix, "%s/", OutDir );
  else prefix[0] = '\0';

  if( R13 < 1.0 )
    snprintf( solv, sizeof solv, "%s%s", "Grad13", Type );
  else
    snprintf( solv, sizeof solv, "%s%s", "Reg13", Type );

  if( !Add ) Convert(t,s2);
  else snprintf( s2, sizeof s2, "All" );

  // common basename incl. output directory, e.g. runs/foo/Reg13_nst1024_Kn0p05_tAll
  if( asprintf(&base, "%s%s_nst%d_Kn%s_t%s", prefix, solv,
      FullGridSize( B, bMax ), Convert(Kn,s1), s2 ) == -1 ) base = NULL;
  if( base == NULL ) { perror("Output: asprintf"); return; }

  if( asprintf(&hdrname, "%s.hdr", base) == -1 ) hdrname = NULL;
  if( asprintf(&datname, "%s.dat", base) == -1 ) datname = NULL;
  if( hdrname == NULL || datname == NULL ) {
    perror("Output: asprintf"); free(base);
    if(hdrname) free(hdrname);
    if(datname) free(datname);
    return;
  };

  // store in the manifest only the .dat file name (strip the directory prefix)
  datbase = strrchr( datname, '/' );
  datbase = datbase ? datbase + 1 : datname;

  if( Add <= 0 ) {                       // fresh files
    Add = 0;
    fh = fopen( hdrname, "w" );
    fd = fopen( datname, "wb" );
    if( fh == NULL || fd == NULL ) {
      perror( fh == NULL ? hdrname : datname );
      if(fh) fclose(fh);
      if(fd) fclose(fd);
      free(base); free(hdrname); free(datname); return;
    };
    WriteManifest( fh, datbase, solv, FullGridSize( B, bMax ), dt, tEnd );
  } else {                               // append a new time instance
    fh = fopen( hdrname, "a" );
    fd = fopen( datname, "ab" );
    if( fh == NULL || fd == NULL ) {
      perror( fh == NULL ? hdrname : datname );
      if(fh) fclose(fh);
      if(fd) fclose(fd);
      free(base); free(hdrname); free(datname); return;
    };
  };

  if( Add == 0 ) printf( "\n  Create File: %s(.hdr/.dat)...", base );
  if( Add > 0 )  printf( "\n  Append (%d) File: %s(.hdr/.dat)...", Add, base );
  fflush(stdout);

  for( nb=0; nb<bMax; nb++ ) {
    lineOut = B[nb].range.rbt.i - B[nb].range.lfb.i + 1;

    // I = points in y-direction (fastest), J = points in x-direction (slowest),
    // identical to the previous Tecplot POINT zone definition.
    ni = B[nb].range.rbt.j - B[nb].range.lfb.j - 2*(MBC-nGhost) + 1;
    nj = B[nb].range.rbt.i - B[nb].range.lfb.i - 2*(MBC-nGhost) + 1;

    // loop start indices (stored in the .hdr so the reader's coordinate
    // reconstruction spans the identical i,j range as this write loop).
    i0 = B[nb].range.lfb.i + MBC - nGhost;
    j0 = B[nb].range.lfb.j + MBC - nGhost;

    nbytes = (long long)ni * nj * NVARS_OUT * (long long)elemSize;

    buf = malloc( (size_t)ni * nj * NVARS_OUT * elemSize );
    if( buf == NULL ) { perror("Output: malloc"); break; }

    // fill the contiguous, variable-fastest buffer in the SAME traversal order
    // as the former POINT loop (i outer over x, j inner over y). x,y are not
    // stored -- only the 10 physical primitives, cast to the output precision.
    p = 0;
    for( i=i0; i<=B[nb].range.rbt.i-MBC+nGhost; i++ ) {
      for( j=j0; j<=B[nb].range.rbt.j-MBC+nGhost; j++ ) {
        dst_idx = (j - B[nb].range.lfb.j)*lineOut + (i - B[nb].range.lfb.i);
        Primitive( B[nb].layer[L_DATA][dst_idx].U, Prim );

        vals[0] = Prim[0];          // rho
        vals[1] = Prim[1];          // vx
        vals[2] = Prim[2];          // vy
        vals[3] = Prim[3]/Prim[0];  // T
        vals[4] = Prim[3];          // p
        vals[5] = Prim[4];          // Pxx
        vals[6] = Prim[5];          // Pxy
        vals[7] = Prim[6];          // Pyy
        vals[8] = Prim[7];          // qx
        vals[9] = Prim[8];          // qy

        v = p*NVARS_OUT;
        if( OutFloat32 )
          for( m=0; m<NVARS_OUT; m++ ) ((float  *)buf)[v+m] = (float)vals[m];
        else
          for( m=0; m<NVARS_OUT; m++ ) ((double *)buf)[v+m] = vals[m];
        p++;
      };
    };

    offset = ftell( fd );
    fwrite( buf, elemSize, (size_t)ni * nj * NVARS_OUT, fd );
    free( buf ); buf = NULL;

    // geometry (BX,BY,DX,DY,LFB_I,LFB_J,ISTART,JSTART,NGHOST) lets the reader
    // rebuild x,y exactly: x = BX + (i-LFB_I)*DX, y = BY + (j-LFB_J)*DY.
    fprintf( fh, "ZONE block=%d I=%d J=%d NVARS=%d "
             "BX=%.16g BY=%.16g DX=%.16g DY=%.16g LFB_I=%d LFB_J=%d "
             "ISTART=%d JSTART=%d NGHOST=%d "
             "SOLUTIONTIME=%.16g STRANDID=%d OFFSET=%ld NBYTES=%lld\n",
             nb+1, ni, nj, NVARS_OUT,
             B[nb].x, B[nb].y, B[nb].dx, B[nb].dy,
             B[nb].range.lfb.i, B[nb].range.lfb.j,
             i0, j0, nGhost,
             t, Add+1, offset, nbytes );
  };

  fclose( fh );
  fclose( fd );
  free( base ); free( hdrname ); free( datname );
  printf( " done\n\n" );

};
