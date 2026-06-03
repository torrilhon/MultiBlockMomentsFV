//*********************************************************************************************
/*                     (struct definitions for 2d R13 solver)                                 */
/*                                                                                            */
/*                                                                                            */
/* (c) Manuel Torrilhon, ETH Zurich, 2010                                                     */
//*********************************************************************************************

#define N 9   // number of variables, i.e., equations

#define NB 10  // maximum number of grid blocks
#define MBC 3  // number of ghost cells


// physical position vector
typedef struct {
  double x, y, z;
} Vector;

// This is the fundamental data structure. It contains the vector of 
// conservative variables in one grid cell.
//   {U[0],U[1],U[2],U[3],U[4],U[5],U[6],U[7],U[8]}
// = { F,  Fx,  Fy,   Fkk,Fxx, Fxy, Fyy, F1x, F1y }
typedef struct {
  double U[N];
} CVector;

// This is an auxiliary data structure that contains the vector of
// primitive variables and its limited derivatives for one grid cell.
typedef struct {
  double Prim[N]; // {rho,vx,vy,p,Pxx,Pxy,Pyy,qx,qy}
  double dPrimX[N]; // {drho,dvx,dvy,dpx,dPxx,dPxy,dPyy,dqx,dqy}
  double dPrimY[N]; // {drho,dvx,dvy,dpy,dPxx,dPxy,dPyy,dqx,dqy}
} AuxData;

// A block represents a rectangular grid with corner (x,y) and spacing (dx,dy)
// and contains three conservative variable vectors for all cells in the form
// of array Data, New, and Tmp. The variable boundary encodes the neighbors
// of the block. (3D is not implemented)
typedef struct {
  CVector *Data, *New, *Tmp;
  AuxData *Aux;
  int block_nr; 
  DB_Range range;
  double x, y, z;
  double dx, dy, dz;
  int boundary[4]; // [0]-left, [1]-right, [2]-bottom, [3]-top, 
                   // extrapolation = 0, wall = -1, synchro = (nr of block)
} DB_Block;

