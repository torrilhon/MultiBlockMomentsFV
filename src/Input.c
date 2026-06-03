//*********************************************************************************************
/*                     (Input-File 2-d R13 equation)                                          */
/*                                                                                            */
/*                                                                                            */
/* (c) Manuel Torrilhon, ETH Zurich, 2010                                                     */
//*********************************************************************************************

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "database_range.h"
#include "defs.h"
#include "tools.h"
#include "Flux.h"


char Type[200];
extern char StartFile[200];



double Kn = 0.001;    // Knudsen number
double F = 0.0;       // volume-force in x-direction
double R13 = 0.0;     // R13 switch


// collision frequency
double CollFreq( double rho, double theta )
{
  double nu0 = 1.0;
  return( nu0*rho*sqrt(theta)/Kn );
};


// local acceleration
void ExternalForce( Vector P, double t, double *g )
{
  double x = P.x, y = P.y;
  
  g[0] = F;
  g[1] = 0.0;
};


//*********************************************************************************************
// provides the initial conditions for the conservative variable vector U at position P.
//*********************************************************************************************
void SetInitial( double *U, Vector P )
{
   double rho, v1, v2, theta, T11, T12, T22, q1, q2;
   double p, vs, E, Prim[N];
   double x = P.x, y = P.y, r = P.x*P.x+P.y*P.y;

  /*** Problem 1  tend=0.3 *****/
  sprintf( Type, "Riemann2D" );

  v1 = 0.0;
  v2 = 0.0;
  T12 = 0.0;
  q1 = 0.0;
  q2 = 0.0;
  theta = 1.0;
  
  if( (sqrt(x*x+y*y) < 0.4) )  {
    rho = 5.0;
  }
  else {
    rho = 1.0;
  };

  T11 = theta;
  T22 = theta;
  /******/

  /*** Problem 2  tend=0.3 *****
  sprintf( Type, "Riemann1D" );

  vx = 0.0;
  vy = 0.0;
  Pxy = 0.0;
  qx = 0.0;
  qy = 0.0;

  if( x < 0.0 )  {
    rho = 5.0;
    vy  = 0.0;
    p   = 5.0;
  }
  else {
    rho = 1.0;
    vy  = 0.0;
    p   = 1.0;
  };

  Pxx = p;
  Pyy = p;
  Pzz = p;
  /******/


  vs  = v1*v1+v2*v2;
  E = 0.5*rho*vs+1.5*rho*theta;

  U[0] = rho;
  U[1] = rho*v1;
  U[2] = rho*v2;
  U[3] = E;
  U[4] = rho*(T11+v1*v1);
  U[5] = rho*(T12+v1*v2);
  U[6] = rho*(T22+v2*v2);
  U[7] = q1+E*v1+rho*(T11*v1+T12*v2);
  U[8] = q2+E*v2+rho*(T12*v1+T22*v2);

}; 


//*********************************************************************************************
// reads block informations from block file 'filename'
// and returns a block list 'B' and output times 'times'
//*********************************************************************************************
int SetBlocks( double *times, DB_Block *B, char *filename )
{
  double xmin, xmax, ymin, ymax, dt, tend;
  int i, n = 0, bMax=0, block_nr, nx, ny, nst = 2*MBC, bound0, bound1, bound2, bound3;
  char tmp[300], *tmp1;
  FILE *fp;
  DB_Range R1;
  
  fp = fopen( filename, "r");
  if( fp == NULL ) { perror(filename); return 0; }

  if( fgets(tmp,300,fp) == NULL ) { fclose(fp); fprintf(stderr,"Unexpected EOF reading times header\n"); return 0; }
  if( fgets(tmp,300,fp) == NULL ) { fclose(fp); fprintf(stderr,"Unexpected EOF reading times header\n"); return 0; }
  
  tmp1 = tmp;
  do times[n++] = strtod(tmp1,&tmp1); while( times[n-1] != 0 );
  if( times[0] > 0.0 ) {
    printf( "   Output times will be: t = " );
    for( i=0; i<n-1; i++ ) printf( "%f, ", times[i] );
  } else {
    dt = -times[0];
    tend = times[1];
    n = (int)((tend+1e-10)/dt);
    for( i=0; i<n; i++ ) times[i] = (i+1)*dt;
    times[n] = 0.0;
    printf( "   %d Output times will be every %f until t_end = %f", n, dt, tend );
  };
    
  if( fgets(tmp,300,fp) == NULL ) { fclose(fp); fprintf(stderr,"Unexpected EOF before block definitions\n"); return 0; }
  if( fgets(tmp,300,fp) == NULL ) { fclose(fp); fprintf(stderr,"Unexpected EOF before block definitions\n"); return 0; }
  printf( "\n\n   Reading block definitions... \n" );
  
  while( fscanf( fp, "%d/%lf,%lf/%lf,%lf/%d,%d/%d,%d,%d,%d/%*[^\n]",
		     &(block_nr), &(xmin), &(xmax), &(ymin), &(ymax), &(nx), &(ny),
         &(bound0), &(bound1), &(bound2), &(bound3) ) != EOF ) {
    printf( "    Block %d, x-range: (%+7.2lf,%+7.2lf), y-range: (%+7.2lf,%+7.2lf), grid %dx%d\n",
		     (block_nr), (xmin), (xmax), (ymin), (ymax), (nx), (ny) );
         
    R1.lfb.i = 1;
    R1.lfb.j = 1;
    R1.rbt.i = nx + nst;
    R1.rbt.j = ny + nst;
  
    bMax++;
    n = block_nr-1;
    B[n].block_nr = block_nr;
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
  int i, j, nb, dst_idx, m, anzMBC = 1;
  double x, y, rho, vx, vy, p, theta, Pxx, Pxy, Pyy, qx, qy, Prim[N];
  char tmp[500];
  Vector P;
  FILE *fp;
    
  sprintf( Type, "FromFile(%s)", filename );

  printf( "\n   Read file: %s...", filename );
  fp = fopen( filename, "r" );
  if( fp == NULL ) { perror(filename); return; }

  if( fgets(tmp,300,fp) == NULL ) { fclose(fp); fprintf(stderr,"Unexpected EOF reading header\n"); return; }
  if( fgets(tmp,300,fp) == NULL ) { fclose(fp); fprintf(stderr,"Unexpected EOF reading header\n"); return; }
  
  for( nb=0; nb<bMax; nb++ ) {
    if( fgets(tmp,300,fp) == NULL ) { fclose(fp); fprintf(stderr,"Unexpected EOF reading block header\n"); return; }

    for( i=B[nb].range.lfb.i+MBC-anzMBC; i<=B[nb].range.rbt.i-MBC+anzMBC; i++ ) {
      for( j=B[nb].range.lfb.j+MBC-anzMBC; j<=B[nb].range.rbt.j-MBC+anzMBC; j++ ) {
        if( fscanf( fp, " %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
                &x, &y, &rho, &vx, &vy, &p, &theta, &Pxx, &Pxy, &Pyy, &qx, &qy ) != 12 ) {
          fclose(fp); fprintf(stderr,"Unexpected data format in %s\n", filename); return; }
        if( fgets(tmp,100,fp) == NULL ) { fclose(fp); fprintf(stderr,"Unexpected EOF after data line\n"); return; }
        printf( "%f %f %f \n", x, y, vx );

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
        Conservativ( Prim, B[nb].Data[dst_idx].U );
      };
    };
  };

  fclose( fp );
  printf( " done\n" );

};
