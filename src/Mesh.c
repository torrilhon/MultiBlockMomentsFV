//*********************************************************************************************
//                     ( grid / restart input: block file + solution reload )
//
//  SetBlocks   -- reads the self-describing block file (key=value header + block table)
//  ParseHeader -- parses the 'key = value' header section (used by SetBlocks)
//  LoadData    -- reloads a solution from a binary (.hdr/.dat) output for restarts
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

double DtOut = 0.0;   // output cadence from the block file's 'dt_out' key; 0 = output
                      // only at the listed 't_out' time(s). Recorded in the header.

int Refine = 1;        // uniform grid refinement factor: multiplies nx,ny of every
                       // block. Set by 'refine' in the block file or the -refine CLI
                       // flag (CLI wins). 1 = grids exactly as written in the file.
int RefineCLI = 0;     // !=0 if -refine was given on the command line; then the
                       // block-file 'refine' value is ignored (CLI overrides it).


//*********************************************************************************************
// parses the global 'key = value' header of the block file. Reads lines from 'fp'
// until the 'Block data' marker (consumed) or EOF. Recognized keys (all optional;
// missing keys keep their built-in defaults; a CLI flag overrides the file value):
//   Kn, R13, F, CFL          -- physics / numerics scalars
//   limiter                  -- none | average | minmod | superbee | harmonic
//   t_out = t1 [t2 ...]      -- requested output time(s); single value or a list
//   dt_out                   -- output cadence (see below)
//   refine                   -- uniform grid refinement factor (default 1)
//
// Output-time semantics (the t=0 initial condition is always written separately):
//   dt_out > 0 : equal intervals of dt_out up to the FIRST t_out value; any
//                further t_out entries are ignored. If t_out/dt_out is non-integer
//                the last interval is shortened so the final output lands
//                exactly on t_out. (-> outputs dt_out, 2dt_out, ..., t_out)
//   dt_out = 0 : no intermediate cadence; output exactly at each listed t_out
//                value. A single t_out therefore gives outputs at t=0 and t=t_out.
//
// Fills 'times' (0-terminated) with the post-IC output times and returns their
// count.
//*********************************************************************************************
static int ParseHeader( FILE *fp, double *times )
{
  char tmp[300], key[64], val[200], *p;
  int i, n, n_times = 0, n_out = 0;
  double dt_out = 0.0, t_out[20] = {0.0};

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
    else if( strcasecmp(key,"dt_out")== 0 ) dt_out = DtOut = atof(val);
    else if( strcasecmp(key,"t_out")== 0 ) {
      char *q = val; double v;
      n_out = 0;
      while( n_out < 19 && (v = strtod(q,&p)) != 0.0 && p != q ) { t_out[n_out++] = v; q = p; }
      t_out[n_out] = 0.0;
    }
    else if( strcasecmp(key,"refine")== 0 ) { if( !RefineCLI ) Refine = atoi(val); }
    else if( strcasecmp(key,"limiter") == 0 ) {
      char name[64]; sscanf(val," %63s",name); Limiter = LimiterCode(name);
    }
    else fprintf( stderr, "   WARNING: unknown parameter '%s' ignored.\n", key );
  };

  // build the post-IC output-time list from dt_out + t_out(list).
  if( dt_out > 0.0 ) {
    // equal cadence up to the FIRST t_out; shorten the last step onto t_out.
    double tfin = t_out[0];
    n = (int)((tfin + 1e-10)/dt_out);          // number of full dt_out steps
    for( i=0; i<n; i++ ) times[i] = (i+1)*dt_out;
    if( n*dt_out < tfin - 1e-10 ) times[n++] = tfin;  // partial final interval
    times[n] = 0.0;
    n_times = n;
    printf( "   Output every dt_out = %g until t_end = %g (%d outputs)\n",
            dt_out, tfin, n_times );
  } else {
    // dt_out = 0: output exactly at each listed t_out value.
    for( i=0; i<n_out; i++ ) times[i] = t_out[i];
    times[n_out] = 0.0;
    n_times = n_out;
    printf( "   Output at t = " );
    for( i=0; i<n_times; i++ ) printf( "%g%s", times[i], i<n_times-1?", ":"\n" );
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
// this loads initial conditions from a single-time binary output (.hdr/.dat)
// pair, matching the layout written by Output(): per block one ZONE record in
// the .hdr (with OFFSET into the .dat) and a variable-fastest payload of the 10
// physical fields (rho vx vy T p Pxx Pxy Pyy qx qy) per point, stored Float32
// (default) or Float64 (PRECISION field). x,y are not stored -- the cells are
// filled in the same i-outer / j-inner traversal order as the writer, so no
// coordinate reconstruction is needed for a restart. 'filename' may be the
// .hdr, the .dat, or the shared basename.
//*********************************************************************************************
void LoadData( DB_Block *B, int bMax, char *filename )
{
  int i, j, nb, dst_idx, nGhost = 1, m, single = 0;
  long zoff;
  double vals[10], Prim[N];
  char base[500], hdrname[520], datname[520], line[600];
  size_t len;
  FILE *fh, *fd;

  // fallback only; the real scenario Type is recovered from the manifest below
  // so a restarted run's output reuses the original scenario name.
  sprintf( Type, "FromFile(%s)", filename );

  // resolve companion paths from .hdr / .dat / bare basename
  strncpy( base, filename, sizeof base - 1 ); base[sizeof base - 1] = '\0';
  len = strlen(base);
  if( len > 4 && ( !strcmp(base+len-4, ".hdr") || !strcmp(base+len-4, ".dat") ) )
    base[len-4] = '\0';
  snprintf( hdrname, sizeof hdrname, "%s.hdr", base );
  snprintf( datname, sizeof datname, "%s.dat", base );

  printf( "\n   Read file: %s(.hdr/.dat)...", base );
  fh = fopen( hdrname, "r" );
  fd = fopen( datname, "rb" );
  if( fh == NULL || fd == NULL ) {
    perror( fh == NULL ? hdrname : datname );
    if(fh) fclose(fh);
    if(fd) fclose(fd);
    return;
  }

  // scan the manifest: pick up PRECISION, then walk it block by block, taking
  // the ZONE records in file order (one per block, matching the write loop).
  nb = 0;
  while( fgets(line, sizeof line, fh) != NULL && nb < bMax ) {
    char *q;
    if( (q = strstr(line, "PRECISION")) != NULL ) {
      single = (strstr(q, "Float32") != NULL);
      continue;
    }
    // recover the scenario name written in the [RUN] section ('Type = ...')
    if( strncmp(line, "Type", 4) == 0 && (q = strchr(line, '=')) != NULL )
      { sscanf( q+1, " %199[^\n]", Type ); continue; }
    if( strncmp(line, "ZONE", 4) != 0 ) continue;

    zoff = 0;
    if( (q = strstr(line, "OFFSET=")) ) zoff = atol(q+7);
    fseek( fd, zoff, SEEK_SET );

    // points in the SAME i-outer (x) / j-inner (y) order the writer used.
    for( i=B[nb].range.lfb.i+MBC-nGhost; i<=B[nb].range.rbt.i-MBC+nGhost; i++ ) {
      for( j=B[nb].range.lfb.j+MBC-nGhost; j<=B[nb].range.rbt.j-MBC+nGhost; j++ ) {
        if( single ) {
          float f[10];
          if( fread(f, sizeof(float), 10, fd) != 10 ) {
            fclose(fh); fclose(fd);
            fprintf(stderr,"Unexpected EOF in %s\n", datname); return; }
          for( m=0; m<10; m++ ) vals[m] = (double)f[m];
        } else {
          if( fread(vals, sizeof(double), 10, fd) != 10 ) {
            fclose(fh); fclose(fd);
            fprintf(stderr,"Unexpected EOF in %s\n", datname); return; }
        }

        // stored physical order: rho vx vy T p Pxx Pxy Pyy qx qy (T redundant).
        Prim[0] = vals[0];   // rho
        Prim[1] = vals[1];   // vx
        Prim[2] = vals[2];   // vy
        Prim[3] = vals[4];   // p   (vals[3] = T = p/rho, recomputed implicitly)
        Prim[4] = vals[5];   // Pxx
        Prim[5] = vals[6];   // Pxy
        Prim[6] = vals[7];   // Pyy
        Prim[7] = vals[8];   // qx
        Prim[8] = vals[9];   // qy
        dst_idx = DBR_Index(i, j, B[nb].range);
        Conservative( Prim, B[nb].layer[L_DATA][dst_idx].U );
      }
    }
    nb++;
  }

  fclose( fh );
  fclose( fd );
  printf( " done\n" );

};
