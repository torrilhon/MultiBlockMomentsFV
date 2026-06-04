/*************************************************************************\
* Definition and operations on database ranges                            *
*                                                                         *
* A database range is a rectangle defined by its two corner points. The   *
* corners are always given by the global indices of the left-bottom and   *
* right-top corner.                                                        *
*                                                                         *
* created:  17.2.99 R.Limacher                                            *
* changed:  moved to 2D only by co-pilot May 2026                         *
\*************************************************************************/

#ifndef _MOT__DBR
#define _MOT__DBR

typedef struct {
  int i, j;
} DB_Point;
/* Data structure addressing a point (cell) in the 2-D grid.
   Indices i, j correspond to the coordinate axes x, y.
*/

typedef struct {
  unsigned int level;
  DB_Point lfb, rbt; /* lfb = left-bottom, rbt = right-top */
} DB_Range;
/* Data structure addressing a rectangular database range.
   level: grid level, 0..DB_MAXLEVEL-1 (see database.c)
   lfb: left-bottom corner of the range
   rbt: right-top corner of the range
   Must satisfy: lfb.i <= rbt.i && lfb.j <= rbt.j
*/


int DBR_Size( DB_Range r);
/* Returns the number of cells in range r.
   r (in): database range
   return = number of cells in r
*/


int DBR_In_Range( int i, int j, DB_Range r);
/* Tests whether cell (i, j) lies within range r.
   i, j (in): indices of the cell
   r (in): database range
   return   0 = cell is not in the range
          <>0 = cell is in the range
*/


int DBR_Index( int i, int j, DB_Range r);
/* Returns the index (number) of cell (i, j) within range r.
   Cell numbering starts at 0 for r.lfb and runs from left to right and
   bottom to top up to DBR_Size(r)-1 for r.rbt.

      -------------
      |   |   |   |
      | 6 | 7 | 8 |
   j  |   |   |   |
   ^  -------------
   |  |   |   |   |
   |  | 3 | 4 | 5 |
   |  |   |   |   |
      -------------
      |   |   |   |
      | 0 | 1 | 2 |
      |   |   |   |
      -------------

         ---> i
 
   i, j (in): indices of the cell
   r (in): database range
   return: index of the cell (if return < 0 or return >= DBR_Size(r), the
           cell is NOT within range r)
*/

#endif
