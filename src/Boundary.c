//*********************************************************************************************
/*                     (boundary routines for R13 esolver in 2d)                              */
/*                                                                                            */
/*                                                                                            */
/* (c) Manuel Torrilhon, ETH Zurich, 2010                                                     */
//*********************************************************************************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "database_range.h"
#include "defs.h"
#include "tools.h"
#include "Input.h"



Vector Get_Point( double i, double j, double k, DB_Block B );

void BoundaryExtrapolate( DB_Block B )
{
  int i, j, n, m;
  int dst_idx, src_idx;
  DB_Range R = B.range;

  for( i=R.lfb.i; i<=R.rbt.i; i++ ) {
    src_idx = DBR_Index(i, R.lfb.j+MBC, B.range); 
    for( n=1; n<=MBC; n++ ) {
      dst_idx = DBR_Index(i, R.lfb.j+MBC-n, B.range);
        for( m=0; m<N; m++ ) {
          B.Data[dst_idx].U[m] = B.Data[src_idx].U[m]; 
          B.Tmp[dst_idx].U[m] = B.Tmp[src_idx].U[m]; 
          B.New[dst_idx].U[m] = B.New[src_idx].U[m]; 
        };
    };
    src_idx = DBR_Index(i, R.rbt.j-MBC, B.range); 
    for( n=1; n<=MBC; n++ ) {
      dst_idx = DBR_Index(i, R.rbt.j-MBC+n, B.range);
        for( m=0; m<N; m++ ) {
          B.Data[dst_idx].U[m] = B.Data[src_idx].U[m]; 
          B.Tmp[dst_idx].U[m] = B.Tmp[src_idx].U[m]; 
          B.New[dst_idx].U[m] = B.New[src_idx].U[m]; 
        };
    };
  };

  for( j=R.lfb.j; j<=R.rbt.j; j++ ) {
    src_idx = DBR_Index(R.lfb.i+MBC, j, B.range); 
    for( n=1; n<=MBC; n++ ) {
      dst_idx = DBR_Index(R.lfb.i+MBC-n, j, B.range);
        for( m=0; m<N; m++ ) {
          B.Data[dst_idx].U[m] = B.Data[src_idx].U[m]; 
          B.Tmp[dst_idx].U[m] = B.Tmp[src_idx].U[m]; 
          B.New[dst_idx].U[m] = B.New[src_idx].U[m]; 
        };
    };
    src_idx = DBR_Index(R.rbt.i-MBC, j, B.range); 
    for( n=1; n<=MBC; n++ ) {
      dst_idx = DBR_Index(R.rbt.i-MBC+n, j, B.range);
        for( m=0; m<N; m++ ) {
          B.Data[dst_idx].U[m] = B.Data[src_idx].U[m]; 
          B.Tmp[dst_idx].U[m] = B.Tmp[src_idx].U[m]; 
          B.New[dst_idx].U[m] = B.New[src_idx].U[m]; 
        };
    };
  };    
};

void BoundarySynchro( DB_Block B1, DB_Block B0, int WallFlag )
{
  int i, j, n, m;
  int dst_idx, src_idx;
  DB_Range R1 = B1.range, R0 = B0.range; 

  if( WallFlag == 2 ) // Bottom of B1(=dst) <- Top of B0(=src)
    for( i=R1.lfb.i; i<=R1.rbt.i; i++ ) 
      for( n=1; n<=MBC; n++ ) {
        src_idx = DBR_Index(i, R0.rbt.j-MBC+1-n, B0.range); 
        dst_idx = DBR_Index(i, R1.lfb.j+MBC-n, B1.range);
        for( m=0; m<N; m++ ) {
          B1.Data[dst_idx].U[m] = B0.Data[src_idx].U[m]; 
          B1.Tmp[dst_idx].U[m] = B0.Tmp[src_idx].U[m]; 
          B1.New[dst_idx].U[m] = B0.New[src_idx].U[m]; 
        };
      };

  if( WallFlag == 3 ) // Top of B1(=dst) <- Bottom of B0(=src)
    for( i=R1.lfb.i; i<=R1.rbt.i; i++ ) 
      for( n=1; n<=MBC; n++ ) {
        src_idx = DBR_Index(i, R0.lfb.j+MBC-1+n, B0.range); 
        dst_idx = DBR_Index(i, R1.rbt.j-MBC+n, B1.range);
        for( m=0; m<N; m++ ) {
          B1.Data[dst_idx].U[m] = B0.Data[src_idx].U[m]; 
          B1.Tmp[dst_idx].U[m] = B0.Tmp[src_idx].U[m]; 
          B1.New[dst_idx].U[m] = B0.New[src_idx].U[m]; 
        };
      };
      
  if( WallFlag == 0 ) // Left of B1(=dst) <- Right of B0(=src)
    for( j=R1.lfb.j; j<=R1.rbt.j; j++ ) 
      for( n=1; n<=MBC; n++ ) {
        src_idx = DBR_Index(R0.rbt.i-MBC+1-n, j, B0.range); 
        dst_idx = DBR_Index(R1.lfb.i+MBC-n, j, B1.range);
        for( m=0; m<N; m++ ) {
          B1.Data[dst_idx].U[m] = B0.Data[src_idx].U[m]; 
          B1.Tmp[dst_idx].U[m] = B0.Tmp[src_idx].U[m]; 
          B1.New[dst_idx].U[m] = B0.New[src_idx].U[m]; 
        };
      };
          
  if( WallFlag == 1 ) // Right of B1(=dst) <- Left of B0(=src)
    for( j=R1.lfb.j; j<=R1.rbt.j; j++ ) 
      for( n=1; n<=MBC; n++ ) {
        src_idx = DBR_Index(R0.lfb.i+MBC-1+n, j, B0.range); 
        dst_idx = DBR_Index(R1.rbt.i-MBC+n, j, B1.range);
        for( m=0; m<N; m++ ) {
          B1.Data[dst_idx].U[m] = B0.Data[src_idx].U[m]; 
          B1.Tmp[dst_idx].U[m] = B0.Tmp[src_idx].U[m]; 
          B1.New[dst_idx].U[m] = B0.New[src_idx].U[m]; 
        };
      };
};

void SetBoundary( DB_Block *B, int bMax, double t )
{
  int n, nb, bSrc;
  
  for( nb=0; nb<bMax; nb++ ) 
    BoundaryExtrapolate( B[nb] );

  for( nb=0; nb<bMax; nb++ ) 
    for( n=0; n<4; n++ ) 
      if( (bSrc = B[nb].boundary[n]) > 0 )
        BoundarySynchro( B[nb], B[bSrc-1], n );
          
};

//*********************************************************************************************
/********* wall routines (only adiabatic!) ***************************************************/
//*********************************************************************************************

// x-Wall, WallFlag pointing into domain
void xFluxWall( double *Prim, double *dPrimX, double *dPrimY, double *F, 
                Vector P, double t, int block, int WallFlag )
{
  #include"primDef.h" // Prim -> rho, v1, v2, theta, T11, T12, T22, q1, q2, vs, E
  double dsig121, dq11, dissip, m111, m121, m221, R11, R12;
  
  dissip = theta/CollFreq(rho,theta)*R13; 
  dsig121 = (5*dPrimY[4]+8*dPrimX[5]-2*dPrimY[6]-3*dPrimY[3])/15;
  dq11    = (5*dPrimX[7]+2*dPrimY[8])/3;

  m121 = -dissip*2*dsig121;
  R11 = -dissip*4.8*dq11;

  v1 = 0.0;
  T12 = 0.0;
  R12 = 0.0;
  q1 = 0.0;
  m111 = 0.0;
  m221 = 0.0;

  // m111,m121,m221,R11,R12 are higher-moment wall-flux terms, kept as scaffolding
  // for a fuller (isothermal / accommodation) wall model but not yet used.
  (void)m111; (void)m121; (void)m221; (void)R11; (void)R12;

  F[0] = rho*v1;
  F[1] = rho*T11+rho*v1*v1;
  F[2] = rho*T12+rho*v1*v2;
  F[3] = q1+1.5*rho*theta*v1+rho*T11*v1+rho*T12*v2+0.5*rho*v1*vs;
  F[4] = 1.2*q1+3*rho*T11*v1+rho*v1*v1*v1;
  F[5] = 0.4*q2+rho*T11*v2+2*rho*T12*v1+rho*v2*v1*v1;
  F[6] = 0.4*q1+rho*T22*v1+2*rho*T12*v2+rho*v1*v2*v2;
  F[7] = 0.5*(7*theta*rho*T11-2*theta*rho*theta + 
         6.4*q1*v1+0.8*q2*v2+4*(rho*T11*v1*v1+rho*T12*v1*v2)+rho*T11*vs+(3*rho*theta+rho*vs)*v1*v1);
  F[8] = 0.5*(7*theta*rho*T12 + 
         2.8*(q1*v2+q2*v1)+2*(rho*T12*vs+(rho*T11+rho*T22)*v1*v2)+rho*T12*vs+(3*rho*theta+rho*vs)*v1*v2);

};

// y-Wall, WallFlag pointing into domain
void yFluxWall( double *Prim, double *dPrimX, double *dPrimY, double *G, 
                Vector P, double t, int block, int WallFlag )
{
  #include"primDef.h" // Prim -> rho, v1, v2, theta, T11, T12, T22, q1, q2, vs, E
  double dsig122, dq22, dissip, m112, m122, m222, R12, R22;
  
  dissip = theta/CollFreq(rho,theta)*R13; 
  dsig122 = (5*dPrimX[6]+8*dPrimY[5]-2*dPrimX[4]-3*dPrimX[3])/15;
  dq22    = (5*dPrimY[8]+2*dPrimX[7])/3;

  m122 = -dissip*2*dsig122;
  R22 = -dissip*4.8*dq22;
  
  v2 = 0.0;
  T12 = 0.0;
  R12 = 0.0;
  q2 = 0.0;
  m222 = 0.0;
  m112 = 0.0;

  // m112,m122,m222,R12,R22 are higher-moment wall-flux terms, kept as scaffolding
  // for a fuller (isothermal / accommodation) wall model but not yet used.
  (void)m112; (void)m122; (void)m222; (void)R12; (void)R22;

  G[0] = rho*v2;
  G[1] = rho*T12+rho*v1*v2;
  G[2] = rho*T22+rho*v2*v2;
  G[3] = q2+1.5*rho*theta*v2+rho*T22*v2+rho*T12*v1+0.5*rho*v2*vs;
  G[4] = 0.4*q2+rho*T11*v2+2*rho*T12*v1+rho*v2*v1*v1;
  G[5] = 0.4*q1+rho*T22*v1+2*rho*T12*v2+rho*v1*v2*v2;
  G[6] = 1.2*q2+3*rho*T22*v2+rho*v2*v2*v2;
  G[7] = 0.5*(7*theta*rho*T12 + 
         2.8*(q1*v2+q2*v1)+2*(rho*T12*vs+(rho*T11+rho*T22)*v1*v2)+rho*T12*vs+(3*rho*theta+rho*vs)*v1*v2);
  G[8] = 0.5*(7*theta*rho*T22-2*theta*rho*theta + 
         6.4*q2*v2+0.8*q1*v1+4*(rho*T22*v2*v2+rho*T12*v1*v2)+rho*T22*vs+(3*rho*theta+rho*vs)*v2*v2);
};






