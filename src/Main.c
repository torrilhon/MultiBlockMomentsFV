//*********************************************************************************************
//                          ( Solver for moment equations )
//
//  build with 'make' (see Makefile at repo root), or directly:
//
//    gcc src/*.c -lm -O3 -o build/moments
//
//  (c) Manuel Torrilhon, ETH Zurich, 2010 - RWTH Aachen 2026
//*********************************************************************************************

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>   // mkdir for -out directory
#include "database_range.h"
#include "defs.h"
#include "tools.h"
#include "Flux.h"
#include "Scenario.h"
#include "Limiter.h"
#include "Mesh.h"
#include "Output.h"


int Steps( DB_Block *B, int bMax, double *t, double tend, double CFLmax, double CFLreg );
double TimeUpdate( DB_Block *B, int bMax, double dt, double t );
double Update( CVector *New, CVector *Data, AuxData *Aux, 
               DB_Block B, CVector *Data0, double dt, double t, int flag );
               
void SetBoundary( DB_Block *B, int bMax, int L, double t );


char StartFile[200] = "";
char OutDir[200]    = "";    // output directory (set via -out), "" = current dir


//*********************************************************************************************
// recursively create directory 'path' (like 'mkdir -p'); existing dirs are fine
//*********************************************************************************************
void MakeDirs( const char *path )
{
  char buf[200]; int i;
  snprintf( buf, sizeof buf, "%s", path );
  for( i=1; buf[i]; i++ )
    if( buf[i]=='/' ) { buf[i]='\0'; mkdir(buf,0755); buf[i]='/'; };
  mkdir( buf, 0755 );
};


//*********************************************************************************************
// main function reads command line parameters, initializes grid and 
// calls 'steps' for evolution until tend with intermediate outputs
// reads blockfile for grid construction and output times.
//*********************************************************************************************
int main(int argc, const char **argv)
{   
  int i, j, k = 1, m, n, nb, bMax, nMax, dst_idx;
  double CFLmax, t, tend;
  double U[N], OutTime[20];
  char BlockFile[200] = "blockfile.txt", limname[64];
  DB_Block B[NB];
  DB_Range R;

  // file/output flags are needed before reading the block file
  if (exists_argument(argc, argv, "-sfile"))
      get_string_argument(argc, argv, "-sfile",StartFile);

  if (exists_argument(argc, argv, "-bfile"))
      get_string_argument(argc, argv, "-bfile",BlockFile);

  if (exists_argument(argc, argv, "-out")) {
      get_string_argument(argc, argv, "-out",OutDir);
      MakeDirs( OutDir );  // create output dir (and parents) if missing
  };

  // -refine must be set BEFORE the grid is built (SetBlocks consumes it); RefineCLI
  // makes it take precedence over the block file's 'refine' value.
  if (exists_argument(argc, argv, "-refine")) {
      Refine = (int)get_float_argument(argc, argv, "-refine");
      if (Refine < 1) Refine = 1;
      RefineCLI = 1;
  };

  // read parameters + grid from the block file (sets Kn, R13, F, CFLreg, Limiter)
  bMax = SetBlocks( OutTime, B, BlockFile );

  // command-line flags override the values read from the block file
  if (exists_argument(argc, argv, "-CFL"))
      CFLreg = (double)get_float_argument(argc, argv, "-CFL");
  if (exists_argument(argc, argv, "-Kn"))
      Kn = (double)get_float_argument(argc, argv, "-Kn");
  if (exists_argument(argc, argv, "-R13"))
      R13 = (double)get_float_argument(argc, argv, "-R13");
  if (exists_argument(argc, argv, "-limiter")) {
      get_string_argument(argc, argv, "-limiter", limname);
      Limiter = LimiterCode(limname);
  };

  printf( "\n Running: Kn = %.4g, R13 = %g, CFL = %.2f, limiter = %s, BlockFile = \"%s\"\n",
          Kn, R13, CFLreg, LimiterName(Limiter), BlockFile );
  printf( "   Output directory: \"%s\"\n", OutDir[0] ? OutDir : "." );
  
  printf( "   Initialize: " );
  for( nb=0; nb<bMax; nb++ ) {
    R = B[nb].range;
    printf( "Block %d (Size %d), ",nb+1,(R.rbt.i-R.lfb.i-2*MBC+1)*(R.rbt.j-R.lfb.j-2*MBC+1) );
    
    for( int l=0; l<N_LAYER; l++ )
      B[nb].layer[l] = malloc( DBR_Size(R)*sizeof(CVector) );
    B[nb].Aux = malloc( DBR_Size(R)*sizeof(AuxData) );
    for( dst_idx=0; dst_idx<DBR_Size(R); dst_idx++ ) 
      for( m=0; m<N; m++ )
        B[nb].layer[L_DATA][dst_idx].U[m] = 0.0;
  }; 
  if( strlen(StartFile) > 0 ) LoadData( B, bMax, StartFile );
  else 
    for( nb=0; nb<bMax; nb++ ) { // initialization without boundary 
      R = B[nb].range;
      int line = R.rbt.i - R.lfb.i + 1;
      for( i=R.lfb.i+MBC; i<R.rbt.i-2; i++ )
        for( j=R.lfb.j+MBC; j<R.rbt.j-2; j++ ) {
          dst_idx = (j - R.lfb.j)*line + (i - R.lfb.i);
          SetInitial( U, GetPoint(i,j,k,B[nb]), 0.0 );
          for( m=0; m<N; m++ )
            B[nb].layer[L_DATA][dst_idx].U[m] = U[m];
        };
    };

  SetBoundary( B, bMax, L_DATA, 0.0 );
  printf( "    done\n" );

  Output( B, bMax, 0, 0.0, -1 );
  
  CFLmax = 1.0;
  
  n = 0;
  nMax = 0;
  t = 0.0;

  while( OutTime[n] > 0 ) {
    tend = OutTime[n];
    nMax += Steps( B, bMax, &t, tend, CFLmax, CFLreg ); 
    Output( B, bMax, nMax, t, nMax );
    n++;
  };

  printf( "\n Done after %d steps.\n\n", nMax  );

  for( nb=0; nb<bMax; nb++ ) {
    for( int l=0; l<N_LAYER; l++ ) free(B[nb].layer[l]);
    free(B[nb].Aux);
  };

  return(0);
}


//*********************************************************************************************
// computes L-infinity difference between New and Data grid functions in 
// the block list B
//*********************************************************************************************
double LinfDifference( DB_Block B )
{
  int idx, m;
  double sum1 = 0.0;

  for( idx=0; idx<DBR_Size(B.range); idx++ ) 
    for( m=0; m<N; m++ )
      sum1 = max(sum1,fabs(B.layer[L_DATA][idx].U[m] - B.layer[L_NEW][idx].U[m]));

  return( sum1 );
};
 
//*********************************************************************************************
// computes the evolution of the grid in block list B from *t to tend with 
// step rejection if CFL > CFLmax and adaptive time steping such that
// CFL = CFLreg. Step information is shown every 'showEvery'
// calls single time step function 'TimeUpdate'
//*********************************************************************************************
int Steps( DB_Block *B, int bMax, double *t, double tend, double CFLmax, double CFLreg )
{
  int dst_idx, nStep = 0, nRej = 0, showEvery = 1, m, nb;
  int gridsize = FullGridSize(B,bMax), minutes = 0, seconds = 0;
  double CFL, dt = 0.0001, perfidx = 0.0, diff = 0.0;
  
  gettimeofday(&time0,NULL);
  while( *t < tend ) {
    nStep++;

    if(!(nStep%showEvery)) {
      printf( " step: %3d...", nStep ); fflush(stdout);
    };
    do {
      CFL = 0.0;
      CFL = TimeUpdate( B, bMax, dt, *t );
    
      if( CFL > CFLmax ) {
        dt = dt*CFLreg/CFL;
        if(!(nStep%showEvery)) 
          printf( "\n   step rejected: CFL = %f (next dt = %f)\n   repeat ...", CFL, dt );
        nRej++;
      } else {
        *t += dt;
        dt = dt*CFLreg/CFL;
        if( *t+dt > tend ) dt = tend-(*t);
      };
    } while( CFL > CFLmax );

    if(!(nStep%showEvery)) {
      if( dt > 0 ) {
        diff = 0.0;
        for( nb=0; nb<bMax; nb++ ) diff = max( diff, LinfDifference( B[nb] ) );
        diff /= dt;
        perfidx = (1.0*GetmSec())/(1.0e-3*gridsize*(nStep+nRej));
        minutes = (int)(1.0e-6*perfidx*gridsize*(tend-(*t))/(60*dt));
        seconds = (int)(1.0e-6*perfidx*gridsize*(tend-(*t))/dt)%60;
      };

      printf( " done.  CFL: %lf, time: %lf, next dt: %lf, diff: %lf, perf = %.2f (ca. %dm%ds)\n", 
	      CFL, *t, dt, diff, perfidx, minutes, seconds );fflush(stdout);
    };

    // copy the whole New layer (interior AND ghost cells) back into Data, so the
    // ghost cells stay valid for the first Update of the next step.
    for( nb=0; nb<bMax; nb++ ) 
      for( dst_idx=0; dst_idx<DBR_Size(B[nb].range); dst_idx++ ) 
        for( m=0; m<N; m++ )
          B[nb].layer[L_DATA][dst_idx].U[m] = B[nb].layer[L_NEW][dst_idx].U[m];
  };

  return( nStep );
};



//*********************************************************************************************
// advances the grid in block list B by one step with dt
// uses SSP2-Runge-Kutta, that is, Heun's method
// calls grid update function 'Update' for every block separately
//*********************************************************************************************
double TimeUpdate( DB_Block *B, int bMax, double dt, double t )
{
  int nb;
  double CFLmax = 0.0, CFL1 = 0.0, CFL2 = 0.0;

  /****** two-dimensional second order Heun ****/
  for( nb=0; nb<bMax; nb++ ) 
    CFL1 = max( CFL1, Update( B[nb].layer[L_TMP], B[nb].layer[L_DATA], B[nb].Aux, B[nb], B[nb].layer[L_DATA], dt, t, 0 ) );
  SetBoundary( B, bMax, L_TMP, t+dt );
  
  for( nb=0; nb<bMax; nb++ ) 
    CFL2 = max( CFL2, Update( B[nb].layer[L_NEW], B[nb].layer[L_TMP], B[nb].Aux, B[nb], B[nb].layer[L_DATA], dt, t+dt, 1 ) );
  SetBoundary( B, bMax, L_NEW, t+dt );
  
  CFLmax = max(CFL1, max(CFL2, CFLmax));
  return( CFLmax );
};


//*********************************************************************************************
// produces an dt-update on a rectangular grid given by B using the data stored in 'Data'
// uses the auxiliar data 'Aux' and calls the limiter 'Limit1', the numerical flux functions
// 'xFlux' and 'yFlux', the production procedures and the boundary fluxes.
// Data0 = NULL : only the residuum is computed from Data and stored in New
// Data0 != NULL : the residuum + Data0 is stored in New
// the Heun flag is used for efficient implementation of Heun's method
// Returns the combined CFL number: dt*lambda_x/dx + dt*lambda_y/dy + CFL_prod
// where CFL_prod = (1 - R(2/3*nu*dt)) / 2 is the relaxation contribution (paper eq.89).
//*********************************************************************************************
double Update( CVector *New, CVector *Data, AuxData *Aux, 
	       DB_Block B, CVector *Data0, double dt, double t, int HeunFlag )
{
  int i, j, k = 1, n, idx, xp_idx, xm_idx, yp_idx, ym_idx;
  double primR[N], primL[N];
  double primT[N], primB[N];
  double U0[N], F[N], G[N];
  double CFL = 0, CFL_xflux = 0, CFL_yflux = 0, CFL_prod = 0;
  double  dtdx = dt/B.dx, dtdy = dt/B.dy;
  DB_Range R = B.range;
  int line = R.rbt.i - R.lfb.i + 1;

  for( i=R.lfb.i+1; i<=R.rbt.i-1; i++ )
    for( j=R.lfb.j+1; j<=R.rbt.j-1; j++ ) {
      idx = (j - R.lfb.j)*line + (i - R.lfb.i);
      ProductionsDirect( Data[idx].U, U0, dt, GetPoint(i,j,k,B), t, 2*HeunFlag+1 );
      Primitive( U0, Aux[idx].Prim ); //initial layer projection
      CFL_prod = max( CFL_prod,
                      ProductionsDirect( Data[idx].U, New[idx].U, dt, GetPoint(i,j,k,B), t, 2*HeunFlag ) );
      if( HeunFlag==1 ) {      // 2nd part Heun-method                
        ProductionsDirect( Data0[idx].U, U0, dt, GetPoint(i,j,k,B), t, -1 );
        for( n=0; n<N; n++ ) 
          New[idx].U[n] += U0[n];
      };
      for( n=0; n<N; n++ ) {
        Aux[idx].dPrimX[n] = 0.0;
        Aux[idx].dPrimY[n] = 0.0;
      };
    };

  for( i=R.lfb.i+1; i<=R.rbt.i-1; i++ )
    for( j=R.lfb.j+1; j<=R.rbt.j-1; j++ ) {
      idx = (j - R.lfb.j)*line + (i - R.lfb.i);
      xp_idx = idx + 1;
      xm_idx = idx - 1;
      yp_idx = idx + line;
      ym_idx = idx - line;
      for( n=0; n<N; n++ ) {
        Aux[idx].dPrimX[n] = Limit1(Aux[idx].Prim[n]-Aux[xm_idx].Prim[n],
                                    Aux[xp_idx].Prim[n]-Aux[idx].Prim[n] )/B.dx;
        Aux[idx].dPrimY[n] = Limit1(Aux[idx].Prim[n]-Aux[ym_idx].Prim[n],
                                    Aux[yp_idx].Prim[n]-Aux[idx].Prim[n] )/B.dy;
      };
    };

  for( i=R.lfb.i+1; i<=R.rbt.i-2; i++ )
    for( j=R.lfb.j+1; j<=R.rbt.j-2; j++ ) {  
      idx = (j - R.lfb.j)*line + (i - R.lfb.i);
      xp_idx = idx + 1;
      yp_idx = idx + line;
                                       // y-Flux:
      for( n=0; n<N; n++ ) {
        primB[n] = Aux[idx].Prim[n] + 0.5*B.dy*Aux[idx].dPrimY[n];
        primT[n] = Aux[yp_idx].Prim[n] - 0.5*B.dy*Aux[yp_idx].dPrimY[n];
      };      
      CFL = yFlux( dtdy, B.dy, primB, Aux[idx].dPrimX, Aux[idx].dPrimY, 
                  primT, Aux[yp_idx].dPrimX, Aux[yp_idx].dPrimY, G );
      CFL_yflux = max(CFL,CFL_yflux);
      
      for( n=0; n<N; n++ ) {	  
        New[idx].U[n]    -= dtdy*G[n];
        New[yp_idx].U[n] += dtdy*G[n];
      };
                                     // x-Flux:
      for( n=0; n<N; n++ ) { 
        primL[n] = Aux[idx].Prim[n] + 0.5*B.dx*Aux[idx].dPrimX[n]; 
        primR[n] = Aux[xp_idx].Prim[n] - 0.5*B.dx*Aux[xp_idx].dPrimX[n];
      }; 
      CFL = xFlux( dtdx, B.dx, primL, Aux[idx].dPrimX, Aux[idx].dPrimY,  
                    primR, Aux[xp_idx].dPrimX, Aux[xp_idx].dPrimY, F ); 
      CFL_xflux = max(CFL,CFL_xflux); 
      
      for( n=0; n<N; n++ ) { 
          New[idx].U[n]    -= dtdx*F[n]; 
          New[xp_idx].U[n] += dtdx*F[n]; 
      }; 
    };

  // Note: wall boundaries are now imposed through reflected ghost cells (see
  // Boundary.c, BC_WALL); the interior fluxes above already use the mirrored
  // ghost data at wall faces, so no separate wall-flux correction is needed.

  if( HeunFlag==1 )      // 2nd part Heun-method
    for( idx=0; idx<DBR_Size(B.range); idx++ ) 
      for( n=0; n<N; n++ )
        New[idx].U[n] = 0.5*New[idx].U[n];
    
  return( dt*CFL_xflux/B.dx+dt*CFL_yflux/B.dy+CFL_prod );
};




