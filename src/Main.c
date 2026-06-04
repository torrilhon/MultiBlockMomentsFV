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
#include "Input.h"


int Steps( DB_Block *B, int bMax, double *t, double tend, double CFLmax, double CFLreg );
double TimeUpdate( DB_Block *B, int bMax, double dt, double t );
double Update( CVector *New, CVector *Data, AuxData *Aux, 
               DB_Block B, CVector *Data0, double dt, double t, int flag );
               
void SetBoundary( DB_Block *B, int bMax, double t );
void Output( DB_Block *B, int bMax, int n, double t, int Add );


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
  
  printf( "   Initialize: " );
  for( nb=0; nb<bMax; nb++ ) {
    R = B[nb].range;
    printf( "Block %d (Size %d), ",nb+1,(R.rbt.i-R.lfb.i-2*MBC+1)*(R.rbt.j-R.lfb.j-2*MBC+1) );
    
    B[nb].Data = malloc( DBR_Size(R)*sizeof(CVector) );
    B[nb].New = malloc( DBR_Size(R)*sizeof(CVector) );
    B[nb].Tmp = malloc( DBR_Size(R)*sizeof(CVector) );
    B[nb].Aux = malloc( DBR_Size(R)*sizeof(AuxData) );
    for( dst_idx=0; dst_idx<DBR_Size(R); dst_idx++ ) 
      for( m=0; m<N; m++ )
        B[nb].Data[dst_idx].U[m] = 0.0;
  }; 
  if( strlen(StartFile) > 0 ) LoadData( B, bMax, StartFile );
  else 
    for( nb=0; nb<bMax; nb++ ) { // Initialisation without boundary 
      R = B[nb].range;
      int line = R.rbt.i - R.lfb.i + 1;
      for( i=R.lfb.i+MBC; i<R.rbt.i-2; i++ )
        for( j=R.lfb.j+MBC; j<R.rbt.j-2; j++ ) {
          dst_idx = (j - R.lfb.j)*line + (i - R.lfb.i);
          SetInitial( U, Get_Point(i,j,k,B[nb]) );
          for( m=0; m<N; m++ )
            B[nb].Data[dst_idx].U[m] = U[m];
        };
    };

  SetBoundary( B, bMax, 0.0 );
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
    free(B[nb].Data);
    free(B[nb].New);
    free(B[nb].Aux);
    free(B[nb].Tmp);
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
      sum1 = max(sum1,fabs(B.Data[idx].U[m] - B.New[idx].U[m]));

  return( sum1 );
};
 
//*********************************************************************************************
// computes the evolution of the grid in block list B from *t to tend with 
// step rejection if CFL > CFLmax and adaptive time steping such that
// CFL = CFLreg. Step information is shown every 'anzShow'
// calls single time step function 'TimeUpdate'
//*********************************************************************************************
int Steps( DB_Block *B, int bMax, double *t, double tend, double CFLmax, double CFLreg )
{
  int dst_idx, anz = 0, anzRej = 0, anzShow = 1, m, nb;
  int gridsize = FullGridSize(B,bMax), minutes = 0, seconds = 0;
  double CFL, dt = 0.0001, perfidx = 0.0, diff = 0.0;
  
  gettimeofday(&time0,NULL);
  while( *t < tend ) {
    anz++;

    if(!(anz%anzShow)) {
      printf( " step: %3d...", anz ); fflush(stdout);
    };
    do {
      CFL = 0.0;
      CFL = TimeUpdate( B, bMax, dt, *t );
    
      if( CFL > CFLmax ) {
        dt = dt*CFLreg/CFL;
        if(!(anz%anzShow)) 
          printf( "\n   step rejected: CFL = %f (next dt = %f)\n   repeat ...", CFL, dt );
        anzRej++;
      } else {
        *t += dt;
        dt = dt*CFLreg/CFL;
        if( *t+dt > tend ) dt = tend-(*t);
      };
    } while( CFL > CFLmax );

    if(!(anz%anzShow)) {
      if( dt > 0 ) {
        diff = 0.0;
        for( nb=0; nb<bMax; nb++ ) diff = max( diff, LinfDifference( B[nb] ) );
        diff /= dt;
        perfidx = (1.0*GetmSec())/(1.0e-3*gridsize*(anz+anzRej));
        minutes = (int)(1.0e-6*perfidx*gridsize*(tend-(*t))/(60*dt));
        seconds = (int)(1.0e-6*perfidx*gridsize*(tend-(*t))/dt)%60;
      };

      printf( " done.  CFL: %lf, time: %lf, next dt: %lf, diff: %lf, perf = %.2f (ca. %dm%ds)\n", 
	      CFL, *t, dt, diff, perfidx, minutes, seconds );fflush(stdout);
    };

    for( nb=0; nb<bMax; nb++ ) 
      for( dst_idx=0; dst_idx<DBR_Size(B[nb].range); dst_idx++ ) 
        for( m=0; m<N; m++ )
          B[nb].Data[dst_idx].U[m] = B[nb].New[dst_idx].U[m];
  };

  return( anz );
};



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

static double eps = 1.0e-10;

// limiter selection (set from the block file key 'limiter' or the -limiter flag;
// see Limiter codes / names in Input.c). Default: harmonic mean.
int Limiter = LIM_HARMONIC;

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
//********************** end (limiter functions) *****************************************


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
    CFL1 = max( CFL1, Update( B[nb].Tmp, B[nb].Data, B[nb].Aux, B[nb], B[nb].Data, dt, t, 0 ) );
  SetBoundary( B, bMax, t+dt );
  
  for( nb=0; nb<bMax; nb++ ) 
    CFL2 = max( CFL2, Update( B[nb].New, B[nb].Tmp, B[nb].Aux, B[nb], B[nb].Data, dt, t+dt, 1 ) );
  SetBoundary( B, bMax, t+dt );
  
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
  double U0[N], F[N], G[N], GW[N], FW[N];
  double CFL = 0, CFL_xflux = 0, CFL_yflux = 0, CFL_prod = 0;
  double  dtdx = dt/B.dx, dtdy = dt/B.dy;
  DB_Range R = B.range;
  int line = R.rbt.i - R.lfb.i + 1;

  for( i=R.lfb.i+1; i<=R.rbt.i-1; i++ )
    for( j=R.lfb.j+1; j<=R.rbt.j-1; j++ ) {
      idx = (j - R.lfb.j)*line + (i - R.lfb.i);
      ProductionsDirect( Data[idx].U, U0, dt, Get_Point(i,j,k,B), t, 2*HeunFlag+1 );
      Primitiv( U0, Aux[idx].Prim ); //initial layer projection
      CFL_prod = max( CFL_prod,
                      ProductionsDirect( Data[idx].U, New[idx].U, dt, Get_Point(i,j,k,B), t, 2*HeunFlag ) );
      if( HeunFlag==1 ) {      // 2nd part Heun-method                
        ProductionsDirect( Data0[idx].U, U0, dt, Get_Point(i,j,k,B), t, -1 );
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
    
  if( B.boundary[0] == -1 ) // wall flux for X-boundary cells (left)
    for( j=R.lfb.j+1; j<=R.rbt.j-1; j++ ) { 
      i = R.lfb.i+MBC; // 'i' is left wall-cell, correcting flux i-1/2  
      idx = (j - R.lfb.j)*line + (i - R.lfb.i);  
      xm_idx = idx - 1;  
      xp_idx = idx + 1;  
      for( n=0; n<N; n++ ) {
        primL[n] = Aux[xm_idx].Prim[n] + 0.5*B.dx*Aux[xm_idx].dPrimX[n];
        primR[n] = Aux[idx].Prim[n] - 0.5*B.dx*Aux[idx].dPrimX[n];        
      };
      xFlux( dtdx, B.dx, primL, Aux[xm_idx].dPrimX, Aux[xm_idx].dPrimY, 
                         primR, Aux[idx].dPrimX, Aux[idx].dPrimY, F );
                         
      for( n=0; n<N; n++ ) {
        Aux[idx].dPrimX[n] = (Aux[xp_idx].Prim[n]-Aux[idx].Prim[n])/B.dx;
        Aux[idx].Prim[n] = 1.5*Aux[idx].Prim[n]-0.5*Aux[xp_idx].Prim[n];
      };                                          
      xFluxWall( Aux[idx].Prim, Aux[idx].dPrimX, Aux[idx].dPrimY, FW, 
                 Get_Point(i-0.5,j,k,B), t, B.block_nr, 1 );

      for( n=0; n<N; n++ )    
        New[idx].U[n] += dtdx*(FW[n]-F[n]);
    };   

  if( B.boundary[1] == -1 ) // wall flux for X-boundary cells (right)
    for( j=R.lfb.j+1; j<=R.rbt.j-1; j++ ) { 
      i = R.rbt.i-MBC; // 'i' is right wall-cell, correcting flux i+1/2  
      idx = (j - R.lfb.j)*line + (i - R.lfb.i);  
      xp_idx = idx + 1;  
      xm_idx = idx - 1;  
      for( n=0; n<N; n++ ) {
        primL[n] = Aux[idx].Prim[n] + 0.5*B.dx*Aux[idx].dPrimX[n];
        primR[n] = Aux[xp_idx].Prim[n] - 0.5*B.dx*Aux[xp_idx].dPrimX[n];        
      };
      xFlux( dtdx, B.dx, primL, Aux[idx].dPrimX, Aux[idx].dPrimY, 
                         primR, Aux[xp_idx].dPrimX, Aux[xp_idx].dPrimY, F );
                         
      for( n=0; n<N; n++ ) {
        Aux[idx].dPrimX[n] = (-Aux[xm_idx].Prim[n]+Aux[idx].Prim[n])/B.dx;
        Aux[idx].Prim[n] = 1.5*Aux[idx].Prim[n]-0.5*Aux[xm_idx].Prim[n];
      };                                          
      xFluxWall( Aux[idx].Prim, Aux[idx].dPrimX, Aux[idx].dPrimY, FW, 
                 Get_Point(i+0.5,j,k,B), t, B.block_nr, -1 );

      for( n=0; n<N; n++ )    
        New[idx].U[n] -= dtdx*(FW[n]-F[n]);
    };   
    
  if( B.boundary[2] == -1 ) // wall flux for Y-boundary cells (bottom)
    for( i=R.lfb.i+1; i<=R.rbt.i-1; i++ ) { 
      j = R.lfb.j+MBC; // 'j' is bottom wall-cell, correcting flux j-1/2  
      idx = (j - R.lfb.j)*line + (i - R.lfb.i);  
      ym_idx = idx - line;  
      yp_idx = idx + line;  
      for( n=0; n<N; n++ ) {
        primB[n] = Aux[ym_idx].Prim[n] + 0.5*B.dy*Aux[ym_idx].dPrimY[n];
        primT[n] = Aux[idx].Prim[n] - 0.5*B.dy*Aux[idx].dPrimY[n];        
      };
      yFlux( dtdy, B.dy, primB, Aux[ym_idx].dPrimX, Aux[ym_idx].dPrimY, 
                         primT, Aux[idx].dPrimX, Aux[idx].dPrimY, G );
                         
      for( n=0; n<N; n++ ) {
        Aux[idx].dPrimY[n] = (Aux[yp_idx].Prim[n]-Aux[idx].Prim[n])/B.dy;
        Aux[idx].Prim[n] = 1.5*Aux[idx].Prim[n]-0.5*Aux[yp_idx].Prim[n];
      };                                          
      yFluxWall( Aux[idx].Prim, Aux[idx].dPrimX, Aux[idx].dPrimY, GW, 
                 Get_Point(i,j-0.5,k,B), t, B.block_nr, 1 );

      for( n=0; n<N; n++ )    
        New[idx].U[n] += dtdy*(GW[n]-G[n]);
    };   

  if( B.boundary[3] == -1 ) // wall flux for Y-boundary cell (top)
    for( i=R.lfb.i+1; i<=R.rbt.i-1; i++ ) { 
      j = R.rbt.j-MBC; // 'j' is top wall-cell, correcting flux j+1/2  
      idx = (j - R.lfb.j)*line + (i - R.lfb.i);  
      yp_idx = idx + line;  
      ym_idx = idx - line;  
      for( n=0; n<N; n++ ) {
        primB[n] = Aux[idx].Prim[n] + 0.5*B.dy*Aux[idx].dPrimY[n];
        primT[n] = Aux[yp_idx].Prim[n] - 0.5*B.dy*Aux[yp_idx].dPrimY[n];        
      };
      yFlux( dtdy, B.dy, primB, Aux[idx].dPrimX, Aux[idx].dPrimY, 
                         primT, Aux[yp_idx].dPrimX, Aux[yp_idx].dPrimY, G );
                         
      for( n=0; n<N; n++ ) {
        Aux[idx].dPrimY[n] = (-Aux[ym_idx].Prim[n]+Aux[idx].Prim[n])/B.dy;
        Aux[idx].Prim[n] = 1.5*Aux[idx].Prim[n]-0.5*Aux[ym_idx].Prim[n];
      };                                          
      yFluxWall( Aux[idx].Prim, Aux[idx].dPrimX, Aux[idx].dPrimY, GW, 
                 Get_Point(i,j+0.5,k,B), t, B.block_nr, -1 );

      for( n=0; n<N; n++ )    
        New[idx].U[n] -= dtdy*(GW[n]-G[n]);
    };   

  if( HeunFlag==1 )      // 2nd part Heun-method
    for( idx=0; idx<DBR_Size(B.range); idx++ ) 
      for( n=0; n<N; n++ )
        New[idx].U[n] = 0.5*New[idx].U[n];
    
  return( dt*CFL_xflux/B.dx+dt*CFL_yflux/B.dy+CFL_prod );
};


//************************************************************
//********** Output Procedures *******************************
//************************************************************

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
  int i, j, k = 1, nb, dst_idx, anzMBC = 1, lineOut;
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
          B[nb].range.rbt.j-B[nb].range.lfb.j-2*(MBC-anzMBC)+1, 
          B[nb].range.rbt.i-B[nb].range.lfb.i-2*(MBC-anzMBC)+1 ); 
    fprintf( fp, ", SOLUTIONTIME=%5.3f, STRANDID=%d\n", t, Add+1 );
    
    for( i=B[nb].range.lfb.i+MBC-anzMBC; i<=B[nb].range.rbt.i-MBC+anzMBC; i++ ) {
      for( j=B[nb].range.lfb.j+MBC-anzMBC; j<=B[nb].range.rbt.j-MBC+anzMBC; j++ ) {
        P = Get_Point(i,j,k,B[nb]);
        dst_idx = (j - B[nb].range.lfb.j)*lineOut + (i - B[nb].range.lfb.i);
        Primitiv( B[nb].Data[dst_idx].U, Prim );

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




