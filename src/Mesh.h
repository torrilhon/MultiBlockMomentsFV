//*********************************************************************************************
//                     ( grid / restart input: block file + solution reload )
//
//  SetBlocks reads the self-describing block file (a 'key = value' header followed
//  by a positional per-block table). LoadData reloads a solution from a Tecplot
//  output file for restarts.
//*********************************************************************************************

extern double CFLreg;   // regulated (target) CFL number

int  SetBlocks( double *times, DB_Block *B, char *filename );
void LoadData( DB_Block *B, int bMax, char *filename );
