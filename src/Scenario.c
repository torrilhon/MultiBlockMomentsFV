//*********************************************************************************************
//                     ( problem definition: physics the user edits )
//
//  This file defines a particular simulation:
//    CollFreq      -- collision frequency (relaxation model)
//    ExternalForce -- body force / acceleration field
//    SetInitial    -- initial conditions for the conservative state
//  To set up a new study you normally only edit this file.
//
//  (c) Manuel Torrilhon, ETH Zurich, 2010 - RWTH Aachen 2026
//*********************************************************************************************

#include <stdio.h>
#include <math.h>
#include "database_range.h"
#include "defs.h"
#include "Scenario.h"


char Type[200];

double Kn  = 0.001;   // Knudsen number
double F   = 0.0;     // volume-force in x-direction
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
  // currently a constant body force in x; position P / time t are available
  // here if a space- or time-dependent force is needed.
  g[0] = F;
  g[1] = 0.0;
};


//*********************************************************************************************
// provides the initial conditions for the conservative variable vector U at position P.
//*********************************************************************************************
void SetInitial( double *U, Vector P )
{
   double rho, v1, v2, theta, T11, T12, T22, q1, q2;
   double vs, E;
   double x = P.x, y = P.y;

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
  //-----

  /*** Problem 2  tend=0.3 -----
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
  -----*/


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
