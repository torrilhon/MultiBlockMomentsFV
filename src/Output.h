//*********************************************************************************************
//                     ( binary output of the solution )
//*********************************************************************************************

extern int OutFloat32;   // 1 = single-precision .dat output (default), 0 = double

char *Convert( double a, char *s );
void  Output( DB_Block *B, int bMax, int n, double t, double dt, double tEnd, int Add );
