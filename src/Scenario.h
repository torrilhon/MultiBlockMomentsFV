//*********************************************************************************************
//                     ( problem definition: physics the user edits )
//
//  Everything in Scenario.c defines a particular simulation: the collision model,
//  body force, and initial state (and, later, time-dependent boundary values).
//  To set up a new study you normally only edit Scenario.c.
//*********************************************************************************************

// physics scalars (set by built-in defaults, the block file, or CLI flags)
extern double Kn;    // Knudsen number
extern double R13;   // R13 switch
extern double F;     // volume force in x-direction
extern char   Type[200];

// specular-wall sign table, indexed [side][variable]; see Scenario.c
extern const double WallSign[4][N];

double CollFreq( double rho, double theta );
void   ExternalForce( Vector P, double t, double *g );
void   SetInitial( double *U, Vector P, double t );
