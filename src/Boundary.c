//*********************************************************************************************
/*                     (boundary routines for the 2d moment solver)                          */
/*                                                                                            */
/*  Ghost-cell boundary treatment. Every edge of every block is filled by one routine,       */
/*  FillEdge(), which handles all four boundary kinds (BCType):                               */
/*                                                                                            */
/*    BC_SYNC      ghost <- interior of the neighbouring block (multi-block coupling)         */
/*    BC_OPEN      ghost <- nearest interior cell        (zero-gradient / extrapolation)      */
/*    BC_WALL      ghost <- mirrored interior cell * sign (specular reflection, WallSign)     */
/*    BC_DIRICHLET ghost <- SetInitial(U, position, t)   (prescribed value, time dependent)   */
/*                                                                                            */
/*  The edge index convention (side) matches DB_Block.boundary[]:                             */
/*    0 = left, 1 = right, 2 = bottom, 3 = top.                                               */
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
#include "Scenario.h"


//---------------------------------------------------------------------------------------------
// copy the full conservative vector (all N components) from cell si of src to cell di of dst
//---------------------------------------------------------------------------------------------
static inline void CopyCell( CVector *dst, int di, const CVector *src, int si )
{
  int m;
  for( m=0; m<N; m++ ) dst[di].U[m] = src[si].U[m];
}

//---------------------------------------------------------------------------------------------
// the opposite edge: left<->right, bottom<->top. Used to locate the facing interior strip of
// a neighbouring block in a synchronisation.
//---------------------------------------------------------------------------------------------
static inline int OppositeSide( int side )
{
  static const int opp[4] = { 1, 0, 3, 2 };
  return opp[side];
}

//---------------------------------------------------------------------------------------------
// (i,j) of the n-th ghost cell (n=1..MBC) on 'side', at tangential index 'ti', in range R.
//---------------------------------------------------------------------------------------------
static inline DB_Point GhostCell( int side, int ti, int n, DB_Range R )
{
  DB_Point p;
  switch( side ) {
    case 0: p.i = R.lfb.i+MBC-n; p.j = ti;            break; // left
    case 1: p.i = R.rbt.i-MBC+n; p.j = ti;            break; // right
    case 2: p.i = ti;            p.j = R.lfb.j+MBC-n; break; // bottom
    case 3: p.i = ti;            p.j = R.rbt.j-MBC+n; break; // top
    default: p.i = ti; p.j = ti;                      break;
  }
  return p;
}

//---------------------------------------------------------------------------------------------
// (i,j) of a source cell on 'side', stepping 's' cells inward from the wall:
//   s = 1   -> the first interior (edge) cell        (used by BC_OPEN)
//   s = n   -> the n-th interior layer               (used by BC_WALL mirror, BC_SYNC)
// Evaluated in range R, which may be the block's own range or a neighbour's.
//---------------------------------------------------------------------------------------------
static inline DB_Point InteriorCell( int side, int ti, int s, DB_Range R )
{
  DB_Point p;
  switch( side ) {
    case 0: p.i = R.lfb.i+MBC-1+s; p.j = ti;              break; // inward from left edge
    case 1: p.i = R.rbt.i-MBC+1-s; p.j = ti;              break; // inward from right edge
    case 2: p.i = ti;              p.j = R.lfb.j+MBC-1+s; break; // inward from bottom edge
    case 3: p.i = ti;              p.j = R.rbt.j-MBC+1-s; break; // inward from top edge
    default: p.i = ti; p.j = ti;                          break;
  }
  return p;
}

//---------------------------------------------------------------------------------------------
// Fill the MBC ghost layers of one 'side' of layer L in block D, according to boundary kind bc.
// For BC_SYNC the source is the neighbouring block S (its facing edge); otherwise S == D.
// 't' is only used by BC_DIRICHLET.
//---------------------------------------------------------------------------------------------
static void FillEdge( DB_Block *D, DB_Block *S, int L, int side, BCType bc, double t )
{
  DB_Range Rd = D->range, Rs = S->range;
  CVector *dst = D->layer[L];
  const double *sgn = (bc==BC_WALL) ? WallSign[side] : NULL;
  int ti, t0, t1, n, di, si, m;
  DB_Point g, s;

  // tangential loop runs along the edge (full range, including corner ghosts)
  if( side < 2 ) { t0 = Rd.lfb.j; t1 = Rd.rbt.j; }   // left/right: vary j
  else           { t0 = Rd.lfb.i; t1 = Rd.rbt.i; }   // bottom/top: vary i

  // conforming-interface check: a sync partner must share the tangential extent
  if( bc == BC_SYNC )
    assert( (side<2) ? (Rd.rbt.j-Rd.lfb.j)==(Rs.rbt.j-Rs.lfb.j)
                     : (Rd.rbt.i-Rd.lfb.i)==(Rs.rbt.i-Rs.lfb.i) );

  for( ti=t0; ti<=t1; ti++ )
    for( n=1; n<=MBC; n++ ) {
      g  = GhostCell( side, ti, n, Rd );
      di = DBR_Index( g.i, g.j, Rd );

      switch( bc ) {
        case BC_SYNC:
          s  = InteriorCell( OppositeSide(side), ti, n, Rs );
          CopyCell( dst, di, S->layer[L], DBR_Index(s.i,s.j,Rs) );
          break;

        case BC_OPEN:
          s  = InteriorCell( side, ti, 1, Rd );           // nearest interior cell
          CopyCell( dst, di, dst, DBR_Index(s.i,s.j,Rd) );
          break;

        case BC_WALL:
          s  = InteriorCell( side, ti, n, Rd );           // mirror the n-th layer
          si = DBR_Index( s.i, s.j, Rd );
          for( m=0; m<N; m++ ) dst[di].U[m] = sgn[m]*dst[si].U[m];
          break;

        case BC_DIRICHLET:
          SetInitial( dst[di].U, GetPoint(g.i, g.j, 1, *D), t );
          break;
      }
    }
}

//---------------------------------------------------------------------------------------------
// map the integer block-file code to a boundary kind
//---------------------------------------------------------------------------------------------
static inline BCType BCFromCode( int code )
{
  if( code >  0 ) return BC_SYNC;
  if( code == 0 ) return BC_OPEN;
  if( code ==-1 ) return BC_WALL;
  return BC_DIRICHLET;   // -2 (and any other negative reserved for future kinds)
}

//---------------------------------------------------------------------------------------------
// Fill the ghost cells of layer L for all blocks, for the four edges each, at time t.
//
// The diagonal MBC x MBC corner ghost blocks are filled WITHOUT any extra corner routine or
// any corner-connectivity data in the block file, by relying on the structure of two sweeps:
//
//   Sweep 1 (all edges, all kinds, full tangential range incl. the corner extension):
//     Every edge ghost strip is written by its own boundary kind. The tangential loop in
//     FillEdge runs the full range, so the corner ghosts are also written -- but as a side
//     effect of whichever edges cover them, in block/edge index order (last writer wins). The
//     corner values from this sweep are therefore order dependent. HOWEVER the *inner*
//     tangential ghost cells (the non-corner part of every strip) are unconditionally correct
//     after sweep 1: open/wall/dirichlet read genuine interior cells of the block, and sync
//     reads genuine interior cells of the neighbour. Interior cells are never modified by
//     SetBoundary, so their values do not depend on processing order.
//
//   Sweep 2 (sync edges only, full tangential range):
//     Re-copy each sync strip from the neighbour. At a corner-tangential position the source
//     index resolves to the neighbour's inner ghost strip -- which sweep 1 already filled
//     correctly -- so the corner ghosts now receive physically valid, order-independent data
//     relayed around the corner. At a re-entrant multi-block notch the vertical-sync and
//     horizontal-sync neighbours fill the corner block differently; that is intentional, since
//     each direction then carries the wall-reflected state appropriate to its own flux.
//
// The only corners left order dependent are convex DOMAIN corners where BOTH adjacent edges
// are non-sync (OO, WW, OW, OD, WD, DD). They keep their finite sweep-1 value; this is a small
// localised inconsistency in the cross-derivative stencil, never a NaN, and is accepted.
//---------------------------------------------------------------------------------------------
void SetBoundary( DB_Block *B, int bMax, int L, double t )
{
  int nb, s, code;
  BCType bc;

  // sweep 1: every edge by its own kind (open/wall/dirichlet from the block itself, sync from
  // the neighbour interior), over the full tangential range including the corner extension.
  for( nb=0; nb<bMax; nb++ )
    for( s=0; s<4; s++ ) {
      bc = BCFromCode( B[nb].boundary[s] );
      if( bc == BC_SYNC ) FillEdge( &B[nb], &B[B[nb].boundary[s]-1], L, s, BC_SYNC, t );
      else                FillEdge( &B[nb], &B[nb],                  L, s, bc,      t );
    }

  // sweep 2: sync edges only. Reads the neighbour's now-correct inner ghost strips into the
  // corner ghosts, giving valid order-independent corner data (no separate corner routine).
  for( nb=0; nb<bMax; nb++ )
    for( s=0; s<4; s++ )
      if( (code = B[nb].boundary[s]) > 0 )
        FillEdge( &B[nb], &B[code-1], L, s, BC_SYNC, t );
}
