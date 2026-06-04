//*********************************************************************************************
//                     ( grid / restart input: block file + solution reload )
//
//  SetBlocks   -- reads the self-describing block file (key=value header + block table)
//  ParseHeader -- parses the 'key = value' header section (used by SetBlocks)
//  LoadData    -- reloads a solution from a Tecplot output file for restarts
//
//  (c) Manuel Torrilhon, ETH Zurich, 2010 - RWTH Aachen 2026
//*********************************************************************************************

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>   // strcasecmp / strncasecmp for the key=value parser
#include "database_range.h"
#include "defs.h"
#include "tools.h"
#include "Flux.h"
#include "Scenario.h"   // Kn, R13, F, Type
#include "Limiter.h"    // Limiter, LimiterCode, LimiterName
#include "Mesh.h"


double CFLreg = 0.95; // regulated CFL number (target)

static int Refine = 1; // uniform grid refinement factor (set by 'refine' in the
                       // block file): multiplies nx,ny of every block. 1 = as-is.


//*********************************************************************************************
// parses the global 'key = value' header of the block file. Reads lines from 'fp'
// until the 'Block data' marker (consumed) or EOF. Recognized keys (all optional;
// missing keys keep their built-in defaults; a CLI flag overrides the file value):
//   Kn, R13, F, CFL          -- physics / numerics scalars
//   limiter                  -- none | average | minmod | superbee | harmonic
//   tend, dt                 -- fixed-cadence output (every dt until tend)
//   times = t1 t2 ...        -- explicit, increasing list of output times
//   refine                   -- uniform grid refinement factor (default 1)
// Fills 'times' (0-terminated) and returns the number of output times.
//*********************************************************************************************
static int ParseHeader( FILE *fp, double *times )
{
  char tmp[300], key[64], val[200], *p;
  int i, n, n_times = 0;
  double dt = 0.0, tend = 0.0;

  printf( "\n   Reading parameters...\n" );
  while( fgets(tmp,300,fp) != NULL ) {
    p = tmp;
    while( *p==' ' || *p=='\t' ) p++;            // skip leading blanks
    if( *p=='\0' || *p=='\n' || *p=='#' ) continue;            // blank / comment
    if( p[0]=='/' && p[1]=='/' ) continue;                     // // comment
    if( strncasecmp(p,"Block data",10) == 0 ) break;           // end of header

    if( sscanf(p," %63[^= \t] = %199[^\n]", key, val ) != 2 ) continue;
    // strip a trailing inline comment from the value
    if( (p = strchr(val,'#')) ) *p = '\0';

    if(      strcasecmp(key,"Kn")  == 0 ) Kn     = atof(val);
    else if( strcasecmp(key,"R13") == 0 ) R13    = atof(val);
    else if( strcasecmp(key,"F")   == 0 ) F      = atof(val);
    else if( strcasecmp(key,"CFL") == 0 ) CFLreg = atof(val);
    else if( strcasecmp(key,"tend")== 0 ) tend   = atof(val);
    else if( strcasecmp(key,"dt")  == 0 ) dt     = atof(val);
    else if( strcasecmp(key,"refine")== 0 ) Refine = atoi(val);
    else if( strcasecmp(key,"limiter") == 0 ) {
      char name[64]; sscanf(val," %63s",name); Limiter = LimiterCode(name);
    }
    else if( strcasecmp(key,"times") == 0 ) {
      char *q = val;
      n_times = 0;
      do times[n_times++] = strtod(q,&q); while( times[n_times-1] != 0.0 );
      n_times--;  // drop the trailing 0 terminator from the count
    }
    else fprintf( stderr, "   WARNING: unknown parameter '%s' ignored.\n", key );
  };

  // build the output-time list: explicit 'times' wins, else 'dt'+'tend' cadence
  if( n_times > 0 ) {
    times[n_times] = 0.0;   // times[] already filled while parsing the 'times' key
    printf( "   Output times: t = " );
    for( i=0; i<n_times; i++ ) printf( "%g%s", times[i], i<n_times-1?", ":"\n" );
  } else {
    n = (int)((tend+1e-10)/dt);
    for( i=0; i<n; i++ ) times[i] = (i+1)*dt;
    times[n] = 0.0;
    n_times = n;
    printf( "   Output every dt = %g until t_end = %g (%d outputs)\n", dt, tend, n );
  };
  if( Refine < 1 ) Refine = 1;
  printf( "   Kn = %.4g, R13 = %g, F = %g, CFL = %.2f, limiter = %s, refine = %d\n",
          Kn, R13, F, CFLreg, LimiterName(Limiter), Refine );

  return( n_times );
}


//*********************************************************************************************
// reads block informations from block file 'filename'
// and returns a block list 'B' and output times 'times'
//*********************************************************************************************
int SetBlocks( double *times, DB_Block *B, char *filename )
{
  double xmin, xmax, ymin, ymax;
  int n = 0, bMax=0, block_id, nx, ny, nst = 2*MBC, bound0, bound1, bound2, bound3;
  char tmp[300];
  FILE *fp;
  DB_Range R1;

  fp = fopen( filename, "r");
  if( fp == NULL ) { perror(filename); return 0; }

  // ----- global header: 'key = value' lines until the 'Block data' marker -----
  ParseHeader( fp, times );

  // ----- per-block table (unchanged positional format) -----
  if( fgets(tmp,300,fp) == NULL ) { fclose(fp); fprintf(stderr,"Unexpected EOF before block definitions\n"); return 0; }
  printf( "\n   Reading block definitions...\n" );
  
  while( fscanf( fp, "%d/%lf,%lf/%lf,%lf/%d,%d/%d,%d,%d,%d/%*[^\n]",
		     &(block_id), &(xmin), &(xmax), &(ymin), &(ymax), &(nx), &(ny),
         &(bound0), &(bound1), &(bound2), &(bound3) ) != EOF ) {
    if( block_id < 1 || block_id > NB ) {
      fprintf( stderr, "   ERROR: block number %d out of range [1,%d=NB]. "
               "Increase NB in defs.h or fix the block file.\n", block_id, NB );
      fclose(fp); return( 0 );
    };
    nx *= Refine;  // uniform refinement: triple/double/... cells per direction
    ny *= Refine;
    printf( "    Block %d, x-range: (%+7.2lf,%+7.2lf), y-range: (%+7.2lf,%+7.2lf), grid %dx%d\n",
		     (block_id), (xmin), (xmax), (ymin), (ymax), (nx), (ny) );
         
    R1.lfb.i = 1;
    R1.lfb.j = 1;
    R1.rbt.i = nx + nst;
    R1.rbt.j = ny + nst;
  
    bMax++;
    n = block_id-1;
    B[n].block_id = block_id;
    B[n].dx = (xmax-xmin)/nx;
    B[n].dy = (ymax-ymin)/ny;
    B[n].x = xmin-(nst/2-0.5)*B[n].dx;
    B[n].y = ymin-(nst/2-0.5)*B[n].dy;
    B[n].range = R1;
    B[n].dz = 1.0;
    B[n].z = 0.0;
    B[n].boundary[0] = bound0;
    B[n].boundary[1] = bound1;
    B[n].boundary[2] = bound2;
    B[n].boundary[3] = bound3;  
  };
  printf( "   ... done.\n\n" );

  return( bMax );
};



//*********************************************************************************************
// this loads initial conditions from a single time output
// file 'filename' in tecplot format according to the output function.
//*********************************************************************************************
void LoadData( DB_Block *B, int bMax, char *filename )
{
  int i, j, nb, dst_idx, nGhost = 1;
  double x, y, rho, vx, vy, p, theta, Pxx, Pxy, Pyy, qx, qy, Prim[N];
  char tmp[500];
  FILE *fp;
    
  sprintf( Type, "FromFile(%s)", filename );

  printf( "\n   Read file: %s...", filename );
  fp = fopen( filename, "r" );
  if( fp == NULL ) { perror(filename); return; }

  if( fgets(tmp,300,fp) == NULL ) { fclose(fp); fprintf(stderr,"Unexpected EOF reading header\n"); return; }
  if( fgets(tmp,300,fp) == NULL ) { fclose(fp); fprintf(stderr,"Unexpected EOF reading header\n"); return; }
  
  for( nb=0; nb<bMax; nb++ ) {
    if( fgets(tmp,300,fp) == NULL ) { fclose(fp); fprintf(stderr,"Unexpected EOF reading block header\n"); return; }

    for( i=B[nb].range.lfb.i+MBC-nGhost; i<=B[nb].range.rbt.i-MBC+nGhost; i++ ) {
      for( j=B[nb].range.lfb.j+MBC-nGhost; j<=B[nb].range.rbt.j-MBC+nGhost; j++ ) {
        if( fscanf( fp, " %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
                &x, &y, &rho, &vx, &vy, &p, &theta, &Pxx, &Pxy, &Pyy, &qx, &qy ) != 12 ) {
          fclose(fp); fprintf(stderr,"Unexpected data format in %s\n", filename); return; }
        if( fgets(tmp,100,fp) == NULL ) { fclose(fp); fprintf(stderr,"Unexpected EOF after data line\n"); return; }

        Prim[0] = rho;
        Prim[1] = vx;
        Prim[2] = vy;   
        Prim[3] = p;   
        Prim[4] = Pxx;
        Prim[5] = Pxy;   
        Prim[6] = Pyy;   
        Prim[7] = qx;   
        Prim[8] = qy;   
        dst_idx = DBR_Index(i, j, B[nb].range);
        Conservative( Prim, B[nb].layer[L_DATA][dst_idx].U );
      };
    };
  };

  fclose( fp );
  printf( " done\n" );

};
