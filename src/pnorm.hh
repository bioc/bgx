#ifndef _PNORM_HH
#define _PNORM_HH

inline double pnorm(double x){ return .5+.5*erf(M_SQRT1_2*x); }

#endif
