//*********************************************************************************************
/*                     (struct definitions for 2d R13 solver)                                 */
/*                                                                                            */
/*                                                                                            */
/* (c) Manuel Torrilhon, ETH Zurich, 2010                                                     */
//*********************************************************************************************

#define N 9   // number of variables, i.e., equations

#define NB 32  // maximum number of grid blocks
#define MBC 3  // number of ghost cells

// (slope-limiter variants LIM_* are defined in Limiter.h)

// the three time-stepping layers held per block; indices into DB_Block.layer[]
enum { L_DATA = 0,    // current solution (start of step)
       L_TMP  = 1,    // intermediate Heun stage
       L_NEW  = 2,    // next solution
       N_LAYER = 3 };

// boundary-condition kinds, derived from the integer code in the block file:
//   code > 0  -> BC_SYNC      (value = id of the neighbouring block)
//   code == 0 -> BC_OPEN      (zero-gradient / extrapolation outflow)
//   code ==-1 -> BC_WALL      (specular reflection with sign flip, see WallSign)
//   code ==-2 -> BC_DIRICHLET (prescribed state from SetInitial(U,P,t))
typedef enum { BC_SYNC, BC_OPEN, BC_WALL, BC_DIRICHLET } BCType;

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
// and contains the conservative variable vectors for all cells in the three
// time-stepping layers layer[L_DATA], layer[L_TMP], layer[L_NEW] (replacing
// the former Data/Tmp/New pointers). The variable boundary encodes the
// neighbours of the block. (3D is not implemented)
typedef struct {
  CVector *layer[N_LAYER];   // [L_DATA] current, [L_TMP] stage, [L_NEW] next
  AuxData *Aux;
  int block_id; 
  DB_Range range;
  double x, y, z;
  double dx, dy, dz;
  int boundary[4]; // [0]-left, [1]-right, [2]-bottom, [3]-top
                   // open = 0, wall = -1, dirichlet = -2, synchro = (nr of block)
} DB_Block;

