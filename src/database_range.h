/*************************************************************************\
* Definition und Operationen von Datenbank-Bereichen                      *
*                                                                         *
* Ein Datenbankbereich ist ein Rechteck, das durch seine beiden Eckpunkte  *
* beschrieben wird. Die Ecke muss immer durch die globalen Indices        *
* der linken, unteren und der rechten, oberen Ecke definiert sein.        *
*                                                                         *
* erstellt:  17.2.99 R.Limacher                                           *
* geaendert: moved to 2D only by co-Pilot May 2026                        *
\*************************************************************************/

#ifndef _MOT__DBR
#define _MOT__DBR

typedef struct {
  int i, j;
} DB_Point;
/* Datenstruktur zur Adressierung eines Punktes ( Zelle) im 2-D-Gitter.
   Indices i, j entsprechen den Koordinatenachsen x, y.
*/

typedef struct {
  unsigned int level;
  DB_Point lfb, rbt; /* lfb = left-bottom, rbt = right-top */
} DB_Range;
/* Datenstruktur zur Adressierung eines rechteckigen Datenbankbereiches.
   level: Stufe des Gitters, 0..DB_MAXLEVEL-1 (siehe database.c)
   lfb: linke, untere Ecke des Bereiches
   rbt: rechte, obere Ecke des Bereiches
   Es muss also gelten: lfb.i <= rbt.i && lfb.j <= rbt.j
*/


int DBR_Size( DB_Range r);
/* Gibt die Anzahl Zellen im Bereich r zurueck.
   r (in): Datenbank-Bereich
   return = Anzahl Zellen in r
*/


int DBR_In_Range( int i, int j, DB_Range r);
/* Prueft, ob die Zelle ( i, j) im Bereich r liegt.
   i, j (in): Indices der Zelle
   r (in): Datenbank-Bereich
   return   0 = Zelle liegt nicht im Bereich
          <>0 = Zelle liegt im Bereich
*/


int DBR_Index( int i, int j, DB_Range r);
/* Gibt den Index (Nummer) der Zelle ( i, j) innerhalb des Bereiches r
   zurueck. Die Nummerierung der Zellen beginnt bei 0 fuer r.lfb und 
   erfolgt von links nach rechts und vorne nach hinten 
   bis DBR_Size( r)-1 fuer r.rbt.

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
 
   i, j (in): Indices der Zelle
   r (in): Datenbank-Bereich
   return: Index der Zelle ( ist return < 0 oder return >= DBR_Size( r), so
           liegt die Zelle NICHT im Bereich r)
*/

#endif
