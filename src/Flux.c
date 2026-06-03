//*********************************************************************************************
/*                     (HLL-Riemann-Solver for G13 equation)                                  */
/*                                                                                            */
/*                                                                                            */
/* (c) Manuel Torrilhon, ETH Zurich, 2010 - RWTH Aachen 2026                                  */
//*********************************************************************************************

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "database_range.h"
#include "defs.h"
#include "tools.h"
#include "Input.h"



//   {U[0],U[1],U[2],U[3],U[4],U[5],U[6],U[7],U[8]}
// = { F,  Fx,  Fy,   Fkk,Fxx, Fxy, Fyy, F1x, F1y }
// Prim[N] -> {rho,vx,vy,theta,Txx,Txy,Tyy,qx,qy}
// dPrimX[N] ->  {drho,dvx,dvy,dthetax,dTxx,dTxy,dTyy,dqx,dqy}
// dPrimY[N] ->  {drho,dvx,dvy,dthetay,dTxx,dTxy,dTyy,dqx,dqy}


//*********************************************************************************************
// computes primitive variables 'Prim' from conservative vector 'U'
//*********************************************************************************************
void Primitiv( double *U, double *Prim )
{
  double rho, v1, v2, theta, T11, T12, T22, q1, q2;
  double vs, E;

  rho = U[0];
  v1 = U[1]/rho;
  v2 = U[2]/rho;
  vs = v1*v1+v2*v2;
  theta  = (2*U[3]/rho-vs)/3;
  E = U[3];
  T11 = U[4]/rho-v1*v1;
  T12 = U[5]/rho-v1*v2;
  T22 = U[6]/rho-v2*v2;
  q1 = U[7]-E*v1-rho*(T11*v1+T12*v2);
  q2 = U[8]-E*v2-rho*(T12*v1+T22*v2);

  Prim[0] = rho;
  Prim[1] = v1;
  Prim[2] = v2;
  Prim[3] = theta;
  Prim[4] = T11;
  Prim[5] = T12;
  Prim[6] = T22;
  Prim[7] = q1;
  Prim[8] = q2;
};

//*********************************************************************************************
// computes conservative variables 'U' from primitive vector 'Prim'
//*********************************************************************************************
void Conservativ( double *Prim, double *U )
{
  #include"primDef.h" // Prim -> rho, v1, v2, theta, T11, T12, T22, q1, q2, vs, E

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

//****************************************************************************
// computes one-dimensional Grad13-Flux in x-direction for 2D variables
//****************************************************************************
void Grad13Flux1D( double *U, double *F )
{
  double rho, v1, v2, T11, T12, T22, T33, q1, q2;
  double vs, theta, sqth, E;
  double M111, M121, M221, M112, M122, R11, R12;

  rho = U[0];
  v1 = U[1]/rho;
  v2 = U[2]/rho;
  vs = v1*v1+v2*v2;
  theta  = (2*U[3]/rho-vs)/3;
  sqth = sqrt(theta);
  E = U[3];
  T11 = U[4]/rho-v1*v1;
  T12 = U[5]/rho-v1*v2;
  T22 = U[6]/rho-v2*v2;
  q1 = U[7]-E*v1-rho*(T11*v1+T12*v2);
  q2 = U[8]-E*v2-rho*(T12*v1+T22*v2);

  q1+E*v1+rho*(T11*v1+T12*v2);
  q2+E*v2+rho*(T12*v1+T22*v2);

  /**** Grad13 closure: ******/
  M111 = 1.2*q1/rho;
  M121 = M112 = 0.4*q2/rho;
  M221 = M122 = 0.4*q1/rho;
  R11 = 5*theta*theta+7*theta*(T11-theta);
  R12 = 7*theta*T12;
  /*****************************************/
  
  // mass balance
  F[0] = rho*v1;
  // momentum balance
  F[1] = rho*(T11+v1*v1);
  F[2] = rho*(T12+v1*v2);
  // energy balance
  F[3] = E*v1+rho*(T11*v1+T12*v2)+q1;
  // pressure tensor balance
  F[4] = rho*(M111+3*T11*v1+v1*v1*v1);  
  F[5] = rho*(M121+T11*v2+2*T12*v1+v1*v2*v1);
  F[6] = rho*(M221+T22*v1+2*T12*v2+v2*v2*v1);
  // heat flux balance
  F[7] = 2*(q1+rho*(T11*v1+T12*v2))*v1+E*v1*v1+rho*(M111*v1+M112*v2)+0.5*rho*(T11*vs+R11);
  F[8] = q2*v1+q1*v2+rho*(v1*v2*(T11+T22)+1.5*T12*vs)+E*v1*v2+rho*(M121*v1+M122*v2)+0.5*rho*R12;        
};


//**********************************************************************
// computes Grad13 flux F for conservative vector U in direction 
// (nx,ny). n is assumed to be normalized. 
//**********************************************************************
void Grad13Flux2D( double *U, double *F, double nx, double ny )
{
  double UR[N], FR[N];
  double tx, ty, s1, s2;

  tx = ny;
  ty = -nx;

  UR[0] = U[0];
  // velocity rotation
  UR[1] = U[1]*nx+U[2]*ny;
  UR[2] = U[1]*tx+U[2]*ty;
  // energy
  UR[3] = U[3];  
  // pressure tensor rotation
  s1 = U[4]*nx+U[5]*ny;
  s2 = U[5]*nx+U[6]*ny;
  UR[4] = s1*nx+s2*ny;
  UR[5] = s1*tx+s2*ty;
  s1 = U[4]*tx+U[5]*ty;
  s2 = U[5]*tx+U[6]*ty;
  UR[6] = s1*tx+s2*ty;
  // heat flux rotation
  UR[7] = U[7]*nx+U[8]*ny;
  UR[8] = U[7]*tx+U[8]*ty;  

  Grad13Flux1D( UR, FR );
  
  F[0] = FR[0];
  // velocity rotation
  F[1] = FR[1]*nx+FR[2]*tx;
  F[2] = FR[1]*ny+FR[2]*ty;
  // energy 
  F[3] = FR[3];
  // pressure tensor rotation
  s1 = FR[4]*nx+FR[5]*tx;
  s2 = FR[5]*nx+FR[6]*tx;
  F[4] = s1*nx+s2*tx;
  F[5] = s1*ny+s2*ty;
  s1 = FR[4]*ny+FR[5]*ty;
  s2 = FR[5]*ny+FR[6]*ty;
  F[6] = s1*ny+s2*ty;
  // heat flux rotation
  F[7] = FR[7]*nx+FR[8]*tx;
  F[8] = FR[7]*ny+FR[8]*ty;
  
};

//*********************************************************************************************
// computes x-Flux from primitive variables and their derivatives
//*********************************************************************************************
static inline double F_Flux( double *Prim, double *dPrimX, double *dPrimY, double *F )
{
  double U[N];
  Conservativ( Prim, U );
  Grad13Flux2D( U, F, 1.0, 0.0 );
  return( 0.0 );
};

//*********************************************************************************************
// computes y-Flux from primitive variables and their derivatives
//*********************************************************************************************
double G_Flux( double *Prim, double *dPrimX, double *dPrimY, double *G )
{
  double U[N];
  Conservativ( Prim, U );
  Grad13Flux2D( U, G, 0.0, 1.0 );
  return( 0.0 );
};

//*********************************************************************************************
// computes damping factors for direct integration of productions terms
//*********************************************************************************************
static inline double Rdamp( double z, int flag )
{
  switch( flag ) {
      case( -1 ): return( (1+0.5*z-0.25*z*z)/(1+z+0.25*z*z) );
      case( 0 ): return( (1+0.5*z-0.25*z*z)/(1+z+0.25*z*z) );
      case( 1 ): return( (1+z)/(1+z+0.25*z*z) );
      case( 2 ): return( (1+0.5*z)/(1+0.5*z) );
      case( 3 ): return( 1.0/(1+0.5*z) );      
  };
  return( 0.0 );
};

//*********************************************************************************************
// computes a relaxation step of the R13-productions from 'U' onto 'New'
//*********************************************************************************************
double ProductionsDirect( double *U, double *New, double dt, Vector X, double t, int flag )
{
  double Prim[N];
  Primitiv( U, Prim );
  #include"primDef.h" // Prim -> rho, v1, v2, theta, T11, T12, T22, q1, q2, vs, E
  double g[2], coll, P[N], dampfac1, dampfac2;

  coll = CollFreq(rho,theta)*dt;

  if( (flag==0) || (flag==2) ) {
    ExternalForce( X, t, g );
    v1 += dt*g[0];
    v2 += dt*g[1];
    vs = v1*v1+v2*v2;
    E = 0.5*rho*vs+1.5*rho*theta;
  };

  dampfac1 = Rdamp( coll, flag );
  dampfac2 = Rdamp( 2.0/3*coll, flag );
  T11 = dampfac1*(T11-theta)+theta;
  T12 = dampfac1*T12;
  T22 = dampfac1*(T22-theta)+theta;
  q1 = dampfac2*q1;
  q2 = dampfac2*q2;

  New[0] = rho;
  New[1] = rho*v1;
  New[2] = rho*v2;
  New[3] = E;
  New[4] = rho*(T11+v1*v1);
  New[5] = rho*(T12+v1*v2);
  New[6] = rho*(T22+v2*v2);
  New[7] = q1+E*v1+rho*(T11*v1+T12*v2);
  New[8] = q2+E*v2+rho*(T12*v1+T22*v2);

  return( 0.5*(1-dampfac2) );
};

//**********************************************************************
// computes Grad13 productions P for conservative vector U.
//**********************************************************************
void Production( double *U, double *P )
{
  double rho, v1, v2, T11, T12, T22, T33, q1, q2;
  double vs, theta, E, nu;

  rho = U[0];
  v1 = U[1]/rho;
  v2 = U[2]/rho;
  vs = v1*v1+v2*v2;
  theta  = (2*U[3]/rho-vs)/3;
  E = U[3];
  T11 = U[4]/rho-v1*v1;
  T12 = U[5]/rho-v1*v2;
  T22 = U[6]/rho-v2*v2;
  q1 = U[7]-E*v1-rho*(T11*v1+T12*v2);
  q2 = U[8]-E*v2-rho*(T12*v1+T22*v2);

  nu = CollFreq(rho, theta);

  // mass/momentum/energy balance
  P[0] = 0.0; 
  P[1] = 0.0;
  P[2] = 0.0;
  P[3] = 0.0;
  // pressure tensor balance
  P[4] = -nu*rho*(T11-theta);  
  P[5] = -nu*rho*T12;
  P[6] = -nu*rho*(T22-theta);
  // heat flux balance
  P[7] = -nu*rho*((T11-theta)*v1+T12*v2)-nu*2*q1/3;
  P[8] = -nu*rho*(T12*v1+(T22-theta)*v2)-nu*2*q2/3;
};

//*********************************************************************************************
// computes a (very coarse) estimation of the maximal sound speed
//*********************************************************************************************
static inline double MaxSound1( double *Prim, double nx, double ny )
{
  #include"primDef.h" // Prim -> rho, v1, v2, theta, T11, T12, T22, q1, q2, vs, E

  return( 2.0*sqrt(5.0/3*theta) );
};

//*********************************************************************************************
// scale computation of diffusive terms for CFL correction
//*********************************************************************************************
static inline double lamDiff( double kaps )
{
  if( kaps > 4.1422533434191076 ) 
    return( (10+4*(kaps+sqrt(6+4*kaps)))/((1+2*kaps)*(1+2*kaps)) );
  if( kaps > 1 ) return( 1.15-0.15*kaps );

  return( 1.0 );
};


//*********************************************************************************************
// computes numerical flux from left and right primitive variables according to HLL
// diffusive terms are simply added inside the flux function
//*********************************************************************************************
static inline double FnumHLLE( double dtdx, double dx, double *PrimL, double *dPrimLX, double *dPrimLY,
				double *PrimR, double *dPrimRX, double *dPrimRY, double *flx, int flag ) 
{
  int n;
  double fR[N], fL[N], UL[N], UR[N];
  double vR, vL, cR = MaxSound1(PrimR,1,0), cL = MaxSound1(PrimL,-1,0);
  double cMax, dScaleL, dScaleR;
  double lamL, lamR, lam, alph;

  if( flag==1 ) {
    vR = PrimR[1]; vL = PrimL[1];
    dScaleL = F_Flux( PrimL, dPrimLX, dPrimLY, fL );
    dScaleR = F_Flux( PrimR, dPrimRX, dPrimRY, fR );
  } else {
    vR = PrimR[2]; vL = PrimL[2];
    dScaleL = G_Flux( PrimL, dPrimLX, dPrimLY, fL );
    dScaleR = G_Flux( PrimR, dPrimRX, dPrimRY, fR );
  };

  cMax = fabs(0.5*(vL+vR))+max(cR,cL);
  cMax = cMax/lamDiff(max(dScaleL,dScaleR)/(cMax*dx));
  lamL = vL-cL;
  lamR = vR+cR;
  lam = (fabs(lamL)-fabs(lamR))/(lamL-lamR);
  alph = (fabs(lamL)-lam*lamL);
  
  Conservativ( PrimL, UL );
  Conservativ( PrimR, UR );
  
  for( n=0; n<N; n++ )
    flx[n] = 0.5*((1+lam)*fL[n]+(1-lam)*fR[n]-alph*(UR[n]-UL[n]));

  return( cMax );
};


//*********************************************************************************************
// numerical flux from left -> right, interface function for 'udate'-procedure
//*********************************************************************************************
double xFlux( double dtdx, double dx, double *primL, double *dPrimLX, double *dPrimLY,
	      double *primR, double *dPrimRX, double *dPrimRY, double *Flux ) 
{
  double cMax;

  cMax = FnumHLLE( dtdx, dx, primL, dPrimLX, dPrimLY, primR, dPrimRX, dPrimRY, Flux, 1 );

  return( cMax );
};

//*********************************************************************************************
// numerical flux from bottom -> top, interface function for 'udate'-procedure
//*********************************************************************************************
double yFlux( double dtdy, double dy, double *primL, double *dPrimLX, double *dPrimLY,
	      double *primR, double *dPrimRX, double *dPrimRY, double *Flux ) 
{
  double cMax;
				      
  cMax = FnumHLLE( dtdy, dy, primL, dPrimLX, dPrimLY, primR, dPrimRX, dPrimRY, Flux, 2 );
  
  return( cMax );
};

