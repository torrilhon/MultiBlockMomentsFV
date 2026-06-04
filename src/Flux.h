/***************************************************************************/
/*                        ( procedures for R13-fluxes    )                 */
/***************************************************************************/


void Primitive( double *U, double *Prim );
void Conservative( double *Prim, double *U );
double xFlux( double dtdx, double dx, double *primL, double *dPrimLX, double *dPrimLY,
	      double *primR, double *dPrimRX, double *dPrimRY, double *Flux ); 
double yFlux( double dtdy, double dy, double *primL, double *dPrimLX, double *dPrimLY,
	      double *primR, double *dPrimRX, double *dPrimRY, double *Flux );
void xFluxWall( double *Prim, double *dPrimX, double *dPrimY, double *F, 
                Vector P, double t, int block, int WallFlag );
void yFluxWall( double *Prim, double *dPrimX, double *dPrimY, double *G, 
                Vector P, double t, int block, int WallFlag );
double ProductionsDirect( double *Prim, double *U, double dt, Vector P, double t, int flag );

