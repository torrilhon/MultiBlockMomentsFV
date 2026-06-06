//*********************************************************************************************
//                     ( grid / restart input: block file + solution reload )
//
//  SetBlocks reads the self-describing block file (a 'key = value' header followed
//  by a positional per-block table). LoadData reloads a solution from a Tecplot
//  output file for restarts.
//*********************************************************************************************

extern double CFLreg;   // regulated (target) CFL number
extern double TOut[21]; // resolved post-IC output times, 0-terminated (no t=0 entry)
extern int    Refine;     // uniform grid refinement factor (block file / -refine)
extern int    RefineCLI;  // set !=0 before SetBlocks to make -refine override the file

int  SetBlocks( DB_Block *B, char *filename );
void LoadData( DB_Block *B, int bMax, char *filename );
