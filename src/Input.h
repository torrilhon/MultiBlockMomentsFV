
extern double Kn; 
extern double R13; 
extern double F;        // volume force in x-direction
extern double CFLreg;   // regulated (target) CFL number
extern int Limiter;     // slope-limiter selection (defined in Main.c)
extern char Type[200];

int LimiterCode( const char *name );
const char *LimiterName( int code );
double CollFreq( double rho, double theta );
void ExternalForce( Vector P, double t, double *g );
void SetInitial( double *U, Vector P );
int SetBlocks( double *times, DB_Block *B, char *filename );
void LoadData( DB_Block *B, int bMax, char *filename );
