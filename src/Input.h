
extern double Kn; 
extern double R13; 
extern char Type[200];

double CollFreq( double rho, double theta );
void ExternalForce( Vector P, double t, double *g );
void SetInitial( double *U, Vector P );
int SetBlocks( double *times, DB_Block *B, char *filename );
void LoadData( DB_Block *B, int bMax, char *filename );
