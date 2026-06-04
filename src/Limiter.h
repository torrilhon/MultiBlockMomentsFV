//*********************************************************************************************
//                     ( slope limiter for second-order reconstruction )
//
//  Self-contained: the limiter variants (LIM_*), the runtime selection (Limiter),
//  the limiter applied to one-sided differences (Limit1), and the name<->code
//  mapping used by the block file / -limiter flag.
//*********************************************************************************************

// slope-limiter variants for the second-order reconstruction
enum LimiterType {
  LIM_NONE     = 0,   // first order (no reconstruction)
  LIM_AVERAGE  = 1,   // unlimited central average
  LIM_MINMOD   = 2,   // minmod
  LIM_SUPERBEE = 3,   // superbee
  LIM_HARMONIC = 4    // harmonic mean (default)
};

extern int Limiter;   // slope-limiter selection (set from block file or -limiter flag)

double      Limit1( double d1, double d2 );           // limited slope from differences (d1,d2)
int         LimiterCode( const char *name );          // name -> LIM_* code
const char *LimiterName( int code );                  // LIM_* code -> name
