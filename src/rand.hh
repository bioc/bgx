/*
 *  This file is part of BGX, the Bayesian Gene eXpression program.
 *  Copyright (c) 2003-2004  Graeme Ambler <graeme@ambler.me.uk>
 *                2006       Ernest Turro <ernest.turro@ic.ac.uk>
 *
 *  BGX is free software; you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License, version 2, as
 *  published by the Free Software Foundation.
 *
 *  BGX is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#ifndef LOCAL_RAND_INCLUDED
#define LOCAL_RAND_INCLUDED

#include <valarray>

#ifdef __INTEL_COMPILER
#include <mathimf.h>
#else
#include <math.h>
#endif

#include <boost/random.hpp>


template<class GenType, typename RealType = double> class Rand;
template<class EngineType, typename RealType = double> class Boost_Wrap;
typedef Rand<Boost_Wrap<boost::lagged_fibonacci4423> > Random;
//typedef Rand<Boost_Wrap<boost::mt19937> > Random;

template<class EngineType, typename RealType>
class Boost_Wrap{
private:
  EngineType Engine;
  boost::uniform_real<RealType> U01;
  boost::variate_generator<EngineType, boost::uniform_real<RealType> > Generator;
public:
  Boost_Wrap(boost::uint32_t seed) : Engine(seed), U01(), Generator(Engine,U01) {}
  
  void seed(boost::uint32_t seed) { Engine.seed(seed); }
  
  RealType operator()() { return Generator(); }
};

template<class GenType, typename RealType>
class Rand{
private:
  GenType unif;
public:
  Rand(boost::uint32_t seed=643774788) : unif(seed), GM_a1(0), GM_a2(0) {}

  void seed(boost::uint32_t seed) { unif.seed(seed); }

  RealType Uniform() { return unif(); }

  RealType Uniform(RealType a, RealType b) { return (b-a)*this->Uniform()+a; }

  // Don't use boost::normal_distribution, as it is slower than mine.
  RealType Normal();

  // Note: in Normal(,) (below) we use std.dev. not variance!!!
  RealType Normal(RealType mu, RealType sigma) { return mu+this->Normal()*sigma; }

  // Note: in TruncNormal(,) (below) we use std.dev. not variance!!!
  RealType TruncNormal(RealType mu, RealType sigma) { 
    // In principle all we want to return is:
    // return (qnorm(this->Uniform(pnorm(-mu/sigma),pnorm((1.0-mu)/sigma))*sigma)+mu); 
    // However, qnorm returns -INF and INF for values at or near 0 and 1.
    // There are also some problems with precision of the qnorm algorithm. 
    // Consequently INF, -INF or negative phis are some times returned.
    // We attempt to deal with this below, but do sometimes still se a negatve
    // phi generated.
    RealType a=pnorm(-mu/sigma);
    RealType b=pnorm((1.0-mu)/sigma);
    RealType u=this->Uniform(a,b);
    // RealType q=qnorm(u);

    if( u < 0.0000000001 ){
      return( 0.0 );
    } else if( (1.0-u) < 0.0000000001 ){
      return( 1.0 );
    } else{
      return( (sigma*qnorm(u))+mu );
    }
  }

  // Don't use boost::exponential_distribution, as it is slower than mine.
  RealType Exponential() { return -log(this->Uniform()); }

  RealType Exponential(RealType lambda) { return this->Exponential()/lambda; }

  RealType Gamma(RealType a);

  RealType Gamma(RealType a, RealType b) { return this->Gamma(a)/b; }

  RealType Cauchy() { return tan(M_PI*(this->Uniform()-.5)); }

  RealType Beta(RealType a, RealType b) { RealType x=this->Gamma(a); RealType y=this->Gamma(b); return x/(x+y); }

  // Dirichlet overwrites the input valarray with the output array.
  void Dirichlet(std::valarray<RealType>& alpha, int length);

  // Ernest's addition
  int Discrete(RealType *cum, int *x, int n) {
    RealType u = this->Uniform();
    if(u <= cum[0]) return x[0];
    for(int i=0; i < n-1; i++) {
      if(u > cum[i] && u <= cum[i+1]) return x[i+1];
    }
    std::cerr << "invalid cumulative distribution. rand.hh:120\n";
    for(int i=0; i < n; i++ ) std::cerr << cum[i] << " ";
    std::cerr << std::endl;
    exit(1);
  }
  // end Ernest's addition
  
  
private:
  RealType GM_a1;
  RealType GM_a2;
  RealType GM_s;
  RealType GM_s2;
  RealType GM_d;
  RealType GM_q0;
  RealType GM_b;
  RealType GM_sigma;
  RealType GM_c;
};

//#include "rand.cc"

/*
 *  This file is part of BGX, the Bayesian Gene eXpression program.
 *  Copyright (c) 2003-2004  Graeme Ambler <graeme@ambler.me.uk>
 *
 *  BGX is free software; you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License, version 2, as
 *  published by the Free Software Foundation.
 *
 *  BGX is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/* Implementation of Ahrens and Dieter (1982) fast algorithm for generating 
   Gamma deviates for a>=1, Ahrens and Dieter (1974) for 0<a<1.  Makes use of 
   Rand::Uniform, Rand::Exponential and Rand::Normal */
template<class GenType, typename RealType>
RealType Rand<GenType, RealType>::Gamma(RealType a)
{
  const RealType Q[9]= {0.0416666664,   // Coefficients of Chebychev polynomial
		      0.0208333723,   // approximation for $q_0$
		      0.0079849875,
		      0.0015746717,
		      -0.0003349403,
		      0.0003340332,
		      0.0006053049,
		      -0.0004701849,
		      0.0001710320};
  const RealType A[9] = {0.333333333,   // Coefficients of Chebychev polynomial
		       -0.249999949,  // approximation for part of $v$
		       0.199999867,
		       -0.166677482,
		       0.142873973,
		       -0.124385581,
		       0.110368310,
		       -0.112750886,
		       0.104089866};
  const RealType E[7] = {1.000000000,   // Coefficients of Chebychev polynomial
		       0.499999994,   // approximation for $\exp(q)-1$
		       0.166666848,
		       0.041664508,
		       0.008345522,
		       0.001353826,
		       0.000247453};

  const RealType FOURRT2=5.65685424949238;  /* $4\sqrt{2}$ */
  const RealType TAU1=-0.71874483771719;
  const RealType EXPM1=0.36787944117144232159;
  
  RealType x;
  RealType t;

  RealType v;
  RealType q;
  RealType e;
  RealType u;

  RealType temp;
  RealType temp2;

  RealType p;

  // Slow but simple algorithm for a<1 
  if (a < 1.0) {
    e = 1.0 + EXPM1 * a;
    for(;;) {
      p = e * this->Uniform();
      if (p >= 1.0) {
	x = -log((e - p) / a);
	if (this->Exponential() >= (1.0 - a) * log(x))
	  return x;
      } 
      else {
	x = exp(log(p) / a);
	if (this->Exponential() >= x)
	  return x;
      }
    }
  }

  // If a==1 then this reduces to an exponential distribution.
  if(a == 1) return this->Exponential();

  // Fast but complicated algorithm for a>1
  if(a != GM_a1){                             // Step 1.
    GM_a1 = a;
    GM_s2 = a-0.5;
    GM_s = sqrt(GM_s2);
    GM_d = FOURRT2-12*GM_s;
  }

  t = this->Normal();                    // Step 2.
  x = GM_s + t/2;
  if(t >= 0) return x*x;

  u = this->Uniform();                   // Step 3.
  if(GM_d*u <= (t*t*t)) return x*x;

  if(a != GM_a2){                             // Step 4.
    GM_a2 = a;
    temp = 1;
    GM_q0 = 0;
    for(int i=0; i<9; ++i){
      temp /= a;
      GM_q0 += Q[i]*temp;
    }
    if(a <= 3.686){
      GM_b = 0.463+GM_s+0.178*GM_s2;
      GM_sigma = 1.235;
      GM_c = 0.195/GM_s - 0.079 + 0.16*GM_s;
    }
    else if(a <= 13.022){
      GM_b = 1.654 + 0.0076*GM_s2;
      GM_sigma = 1.68/GM_s + 0.275;
      GM_c = 0.062/GM_s + 0.024;
    }
    else{
      GM_b = 1.77;
      GM_sigma = 0.75;
      GM_c = 0.1515/GM_s;
    }
  }

  if(x > 0){                               // Step 5.

    v = t/(GM_s+GM_s);                           // Step 6.
    if((v>0.25) || (v<0.25)) q = GM_q0-GM_s*t+t*t/4+(GM_s2+GM_s2)*log(1+v);
    else{
      temp=1;
      temp2=0;
      for(int i=0; i<9; ++i){
	temp *= v;
	temp2 += A[i]*temp;
      }
      q = GM_q0 + t*t*temp2/2;
    }

    if(log(1-u) <= q) return x*x;          // Step 7.
  }

  do {
    do {
      e = this->Exponential();           // Step 8.
      u = this->Uniform();
      u = u + u - 1;
      if(u > 0) t = GM_b + e*GM_sigma;
      else t = GM_b - e*GM_sigma;
      
    } while (t <= TAU1);                   // Step 9.

    v = t/(GM_s+GM_s);                           // Step 10.
    if((v>0.25) || (v<0.25)) q = GM_q0-GM_s*t+t*t/4+(GM_s2+GM_s2)*log(1+v);
    else{
      temp=1;
      temp2=0;
      for(int i=0; i<9; ++i){
	temp *= v;
	temp2 += A[i]*temp;
      }
      q = GM_q0 + t*t*temp2/2;
    }

    temp2 = 0;                             // Step 11.
    if(q > 0){
      if(q <= 0.5){
	temp=1;
	temp2=0;
	for(int i=0; i<7; ++i){
	  temp *= q;
	  temp2 += E[i]*temp;
	}
	temp2 *= exp(e-t*t/2);
      }
      else
	temp2 = (exp(q)-1)*exp(e-t*t/2);
    }
  } while((q <= 0) || ((u>=0)&&(GM_c*u>temp2)) || ((u<0)&&(GM_c*u*-1>temp2)));

  x=GM_s+t/2;                                 // Step 12.
  return x*x;
}

template<class GenType, typename RealType>
void Rand<GenType, RealType>::Dirichlet(std::valarray<RealType>& a, int length)
{
  RealType sum=0;

  for(int i=0; i<length; i++){
    a[i]=this->Gamma(a[i]);
    sum += a[i];
  }
  for(int i=0; i<length; i++)
    a[i]=a[i]/sum;
}

// Generates Normal() deviates using the Trapezoidal method of Ahrens and 
// Dieter (1972). Refer to the paper for the meaning of constants A--S.
// Very fast, and the samples are nicely Normal.
template<class GenType, typename RealType>
inline RealType Rand<GenType, RealType>::Normal()
{
  const RealType A=0.919544405706926; const RealType B=2.40375765693742;
  const RealType C=0.825339282536923; const RealType D=2.11402808333742;
  const RealType E=0.965487131213858; const RealType F=4.46911473713927;
  const RealType G=0.949990708733028; const RealType H=1.84039874739771;
  const RealType I=0.273629335939706; const RealType J=0.398942280401433;
  const RealType K=0.443299125820220; const RealType L=0.209694057195486;
  const RealType M=0.042702581590795; const RealType N=0.925852333707704;
  const RealType O=0.289729573600000; const RealType P=1.55066917379771;
  const RealType Q=0.015974522655238; const RealType R=0.382544556042518;
  const RealType S=0.016397724358915;

  RealType u=this->Uniform(); RealType u0=this->Uniform();
  if(u<A) return B*(u0+u*C)-D;
  else{
    RealType y;
    RealType u1; RealType u2;
    if(u>=E){
      do{
	u1=this->Uniform(); u2=this->Uniform();
	y=sqrt(F-2*log(u1));
      } while(y*u2 > D);
      if(u0<.5) return y; return -y;
    }
    else if(u>=G){
      do{
	u1=this->Uniform(); u2=this->Uniform();
	y=H+u1*I;
      } while(J*exp(-y*y*.5) - K + y*L < u2*M);
      if(u0<.5) return y; return -y;
    }
    else if(u>=N){
      do{
	u1=this->Uniform(); u2=this->Uniform();
	y=O+u1*P;
      } while(J*exp(-y*y*.5) - K + y*L < u2*Q);
      if(u0<.5) return y; return -y;
    }
    else{
      do{
	u1=this->Uniform(); u2=this->Uniform();
	y=u1*O;
      } while(J*exp(-y*y*.5) - R < u2*S);
    if(u0<.5) return y; return -y;
    }
  }
}

#endif // LOCAL_RAND_INCLUDED
