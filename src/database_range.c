/*************************************************************************\
* Definition und Operationen von Datenbank-Bereichen                      *
*                                                                         *
* erstellt:  17.2.99 R.Limacher                                           *
* geaendert:                                                              *
\*************************************************************************/

#include "database_range.h"

int DBR_Size( DB_Range r) {

  return ( (r.rbt.i - r.lfb.i + 1) *
           (r.rbt.j - r.lfb.j + 1) );
}

int DBR_In_Range( int i, int j, DB_Range r) {

  return ( (i >= r.lfb.i) && (i <= r.rbt.i) &&
           (j >= r.lfb.j) && (j <= r.rbt.j) );
}

int DBR_Index( int i, int j, DB_Range r) {

  int line = r.rbt.i - r.lfb.i + 1;

  return ( (j - r.lfb.j)*line + (i - r.lfb.i) );
}

