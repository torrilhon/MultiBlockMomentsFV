/***************************************************************************/
/*                        ( procedures for R13-fluxes    )                 */
/***************************************************************************/


void Primitive( double *U, double *Prim );
void Conservative( double *Prim, double *U );
double xFlux( double dtdx, double dx, double *primL, double *dPrimLX, double *dPrimLY,
	      double *primR, double *dPrimRX, double *dPrimRY, double *Flux ); 
double yFlux( double dtdy, double dy, double *primL, double *dPrimLX, double *dPrimLY,
	      double *primR, double *dPrimRX, double *dPrimRY, double *Flux );
// (wall fluxes xFluxWall/yFluxWall removed: walls are now imposed via reflected
//  ghost cells in Boundary.c, BC_WALL)
double ProductionsDirect( double *Prim, double *U, double dt, Vector P, double t, int flag );

