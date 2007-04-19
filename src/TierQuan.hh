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

#include <math.h>
#include <vector>
#include <algorithm>
using namespace std;

// p is a number BETWEEN 0 AND 1.
class Quantiles_T{
private:
  vector<double> p;
  int n;
  int init_size;
  vector<double> init;
  vector<double> quantiles;

  double h;
  vector<double> f;
  vector<double> d;
  double a;
  double d_0;

  double quan(double p_) // Find the sample quantiles during initialization
  {
    double location=(n-1)*p_;
    int loc=(int)location;
    nth_element(init.begin(), init.begin()+loc+1, init.begin()+n);
    location-=loc;
    return location*(init[loc+1]-init[loc])+init[loc];
  }

  void update_h() { h=1/sqrt(n-init_size+1.0); }

  void update_f(double x, int i)
  {
    update_h();
    int I=0;
    if( fabs(quantiles[i]-x) <= h ) I=1;
    f[i] = ( (n-init_size)*f[i] + I/(2*h) ) / (n-init_size+1);
  }

  void update_d(double x, int i)
  {
    d[i] = fmin(1/f[i],d_0*pow((double)(n+1),a));
  }

  void update_quantile(double x, int i)
  {
    update_f(x,i);
    int Z=0;
    if( x <= quantiles[i] ) Z=1;
    quantiles[i] -= d[i]*(Z-p[i])/(n-init_size+1);
    update_d(x,i);
  }

public:
  Quantiles_T() {}

  Quantiles_T(vector<double> p_, int init_size_=100)
    : p(p_), n(0), init_size(init_size_), init(init_size_), 
      quantiles(p_.size()), h(1), f(p_.size()), d(p_.size()), a(0.4)
  {
    for(unsigned i=0; i<p_.size(); ++i) { f[i]=0.0000001; }
  }

   void Initialize(vector<double> p_, int init_size_=100)
  {
    p=p_; n=0; init_size=init_size_; init.resize(init_size_); 
    quantiles.resize(p_.size()); h=1; f.resize(p_.size()); 
    d.resize(p_.size()); a=0.4;
    for(unsigned i=0; i<p_.size(); ++i) { f[i]=0.0000001; }
  }

  void operator()(double x)
  {
    if(n>init_size)
      {
	++n;
	for(unsigned i=0; i<p.size(); ++i)
	  update_quantile(x,i);
      }
    else if(n<init_size)
      {
	init[n++]=x;
      }
    else // if (n==init_size)
      {
	if(n>1){
	  for(unsigned i=0; i<p.size(); ++i)
	    quantiles[i]=quan(p[i]);
	  d_0=quan(0.75)-quan(0.25);
	  ++n;
	  for(unsigned i=0; i<p.size(); ++i)
	    update_quantile(x,i);
	}
	else if(n==1){
	  for(unsigned i=0; i<p.size(); ++i)
	    quantiles[i] = p[i]*fabs(x-init[0])+fmin(init[0],x);
	  d_0=0.5*fabs(x-init[0]);
	  ++n;
	}
	else{ //n==0
	  for(unsigned i=0; i<p.size(); ++i)
	    quantiles[i]=x;
	  d_0=x;
	  ++n;
	}
      }
  }

  vector<double> Output()
  {
    if(n==0)
       for(unsigned i=0; i<p.size(); ++i)
	 quantiles[i]=0;
    else if(n<=init_size)
      for(unsigned i=0; i<p.size(); ++i)
	quantiles[i]=quan(p[i]);
    return quantiles;
  }
};

