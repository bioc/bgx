/*
 *  This file is part of BGX, the Bayesian Gene eXpression program.
 *  Copyright (c) 2003-2004  Graeme Ambler <graeme@ambler.me.uk>
 *                2005-2006  Ernest Turro <ernest.turro@ic.ac.uk>
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

#ifndef _BGX_UPDATES_HH
#define _BGX_UPDATES_HH

#include "pnorm.hh"
#include "qnorm.h"
#include "rwm.hh"
#include <valarray>

using namespace std; 
typedef valarray<double> varray;
typedef valarray<valarray<double> > varray2d;

#define SAFE_EXP(a) exp(max(-500.0,min(0.0,a)))

class S_T{
private:
  varray2d& PM;
  varray2d& MM;
  varray2d& H;
  double &phi;
  varray2d& mu;
  varray2d& taug;
  varray& tau;
  varray& beta;
  int* probesets;
  int* samplesets;
  int* categories;
  int g;
  int p;
  int c;
  int r;
public:
  S_T(varray2d& PM_, varray2d& MM_, varray2d& H_, double &phi_, varray2d& mu_, 
      varray2d& taug_, varray& tau_, varray& beta_, int* probesets_,
      int* samplesets_, int * categories_)
    : PM(PM_), MM(MM_), H(H_), phi(phi_), mu(mu_), taug(taug_), 
      tau(tau_), beta(beta_), probesets(probesets_), samplesets(samplesets_), categories(categories_) {}

  double operator()(varray2d& S, double S2, int j, int i)
  {
    if(i==0){
      g=p=0;
      if(j==0) c=r=0;
      else if(++r==samplesets[c]) ++c,r=0;
    }
    else if(++p==probesets[g]) ++g,p=0;
    if(S2>0)
      {
	double x=log(S[j][i]+1);
	double x2=log(S2+1);
	return SAFE_EXP(taug[c][g]*((x-mu[c][g])*(x-mu[c][g])-(x2-mu[c][g])*(x2-mu[c][g]))*.5
			+ (x-x2)
			+ tau[j]*((PM[j][i]-S[j][i]-H[j][i]-beta[j])*(PM[j][i]-S[j][i]-H[j][i]-beta[j]) 
			       + (MM[j][i]-phi*S[j][i]-H[j][i]-beta[j])*(MM[j][i]-phi*S[j][i]-H[j][i]-beta[j])
			       - (PM[j][i]-S2-H[j][i]-beta[j])*(PM[j][i]-S2-H[j][i]-beta[j]) 
			       - (MM[j][i]-phi*S2-H[j][i]-beta[j])*(MM[j][i]-phi*S2-H[j][i]-beta[j]))*.5 );
      }
    else return 0;
  }
};

class H_T{
private:
  varray2d& PM;
  varray2d& MM;
  varray2d& S;
  double &phi;
  varray2d& lambdas;
  //varray2d& etas;
  varray& eta;
  varray& tau;
  varray& beta;
  int * categories; 
public:
  H_T(varray2d& PM_, varray2d& MM_, varray2d& S_, double &phi_, varray2d& lambdas_,
      varray& eta_, varray& tau_, varray& beta_, int *categories_)
    : PM(PM_), MM(MM_), S(S_), phi(phi_), lambdas(lambdas_), eta(eta_), 
      tau(tau_), beta(beta_), categories(categories_) {}
  
  double operator()(varray2d& H, double H2, int j, int i)
  { 
    if(H2>0)
      {
  double y=log(H[j][i]+1);
  double y2=log(H2+1);
  return SAFE_EXP(eta[j]*( (y-lambdas[j][categories[i]])*(y-lambdas[j][categories[i]]) -
             (y2-lambdas[j][categories[i]])*(y2-lambdas[j][categories[i]]) )*.5
      + (y-y2)
      + tau[j]*((PM[j][i]-S[j][i]-H[j][i]-beta[j])*(PM[j][i]-S[j][i]-H[j][i]-beta[j])
             + (MM[j][i]-phi*S[j][i]-H[j][i]-beta[j])*(MM[j][i]-phi*S[j][i]-H[j][i]-beta[j])
             - (PM[j][i]-S[j][i]-H2-beta[j])*(PM[j][i]-S[j][i]-H2-beta[j]) 
             - (MM[j][i]-phi*S[j][i]-H2-beta[j])*(MM[j][i]-phi*S[j][i]-H2-beta[j]))*.5 );
      } 
    else return 0;
  }
};

class Mu_T{
private:
  varray2d& S;
  varray2d& taug;
  int* probesets;
  int* samplesets;
  int offset;
  int offset2;
public:
  Mu_T(varray2d& S_, varray2d& taug_, int* probesets_, int* samplesets_)
    : S(S_), taug(taug_), probesets(probesets_), samplesets(samplesets_) {}

  double operator()(varray2d& mu, double mu2, int k, int i)
  {
    if(i==0){
      offset=0;
      if(k==0) offset2=0;
      else offset2 += samplesets[k-1];
    }
    else offset += probesets[i-1];
    if(mu2>0 && mu2<15)
      {
	double sum=0;
	for(int l=0; l<samplesets[k]; ++l)
	  {
	    for(int j=0; j<probesets[i]; ++j)
	      {
		double x=log(S[offset2+l][offset+j]+1);
		sum += (x-mu[k][i])*(x-mu[k][i]) - (x-mu2)*(x-mu2);
	      }
	  }
	return SAFE_EXP(taug[k][i]*( sum )*.5 + 
			samplesets[k]*probesets[i]*log(pnorm(mu[k][i]*sqrt(taug[k][i]))/
						       pnorm(mu2*sqrt(taug[k][i]))));
      }
    else return 0;
  }
};

class Sigma_T{
private:
  varray2d& S;
  varray2d& mu;
  varray& a;
  varray& b;
  int* probesets;
  int* samplesets;
  int offset;
  int offset2;
public:
  Sigma_T(varray2d& S_, varray2d& mu_, varray& a_, varray& b_, int* probesets_, int* samplesets_)
    : S(S_), mu(mu_), a(a_), b(b_), probesets(probesets_), samplesets(samplesets_) {}

  double operator()(varray2d& taug, double taug2, int k, int i)
  {
    if(i==0){
      offset=0;
      if(k==0) offset2=0;
      else offset2 += samplesets[k-1];
    }
    else offset += probesets[i-1];
    if(taug2>0)
      {
	double sum=0;
	for(int l=0; l<samplesets[k]; ++l)
	  {
	    for(int j=0; j<probesets[i]; ++j)
	      {
		double x=log(S[offset2+l][offset+j]+1);
		sum += (x-mu[k][i])*(x-mu[k][i]);
	      }
	  }
	return 
	  SAFE_EXP((taug[k][i]-taug2)*( sum )*.5
		   + (samplesets[k]*probesets[i]*.5-1)*log(taug2/taug[k][i])
		   + samplesets[k]*probesets[i]*log(pnorm(mu[k][i]*sqrt(taug[k][i]))/
			       pnorm(mu[k][i]*sqrt(taug2)))
		   + b[k]*((log(taug[k][i])-a[k])*(log(taug[k][i])-a[k])
			   -(log(taug2)-a[k])*(log(taug2)-a[k]))*.5);
      }
    else return 0;
  }
};


class Lambda_T{
private:
  varray2d& H;
  //varray2d& etas;
  varray& eta;
  double p_mu;
  double p_tau;
  int probes;
  vector<vector<int> > &cat_indices;
public:
  Lambda_T(varray2d& H_, varray& eta_, double p_mu_, double p_tau_, vector<vector<int> > & cat_indices_)
    : H(H_), eta(eta_), p_mu(p_mu_), p_tau(p_tau_), probes(H[0].size()),cat_indices(cat_indices_){}
  
  double operator()(varray2d& lambda, double lambda2, int j, int cat) //j == samples
  {
    double sum=0;
    for(unsigned int i=0; i < cat_indices[cat].size(); ++i)
      { 
  double yi=log(H[j][cat_indices[cat][i]]+1);
  sum += (yi-lambda[j][cat])*(yi-lambda[j][cat]) - (yi-lambda2)*(yi-lambda2);
      }
    
    return SAFE_EXP(eta[j]*sum*.5 + /*probes*/cat_indices[cat].size()*log(pnorm(lambda[j][cat]*sqrt(eta[j]))/
                 pnorm(lambda2*sqrt(eta[j])))
        + p_tau*( (lambda[j][cat]-p_mu)*(lambda[j][cat]-p_mu)
            - (lambda2-p_mu)*(lambda2-p_mu) )*.5 );
  }
};


class Eta_T{
private:
  varray2d& H;
  varray2d& lambdas;
  double alpha;
  double beta;
  int probes;
  vector<vector<int> > &cat_indices;
  int *categories;
  int noOfCategories;
public:
  Eta_T(varray2d& H_, varray2d& lambdas_, double alpha_, double beta_, vector<vector<int> > &cat_indices_,int *categories_)
    : H(H_), lambdas(lambdas_), alpha(alpha_), beta(beta_), probes(H[0].size()), cat_indices(cat_indices_), 
    categories(categories_),noOfCategories(lambdas[0].size()){}

  
double operator()(varray & eta, double eta2, int j/*, int cat*/)
  {
    double sum=0.0;
    for(int i=0;i < probes; ++i) {
      sum += (log(H[j][i]+1)-lambdas[j][categories[i]])*(log(H[j][i]+1)-lambdas[j][categories[i]]);
    }
    double sum2=0.0;
    for(int i=0; i < noOfCategories; i++) {
      sum2 += cat_indices[i].size()*log(pnorm(lambdas[j][i]*sqrt(eta[j]))/pnorm(lambdas[j][i]*sqrt(eta2)));
    }

    return SAFE_EXP((eta[j]-eta2)*(.5*sum+beta)
        + sum2
        + (alpha+.5*probes-1)*log(eta2/eta[j]));
  }
};

class Kappa_T{
private:
  varray2d& lambdas;
  double alpha;
  double beta;
  int probes;
  vector<vector<int> > &cat_indices;
  int noOfCategories;
  int *categories;
public:
  Kappa_T(varray2d& lambdas_, double alpha_, double beta_, int probes_, vector<vector<int> > &cat_indices_,int *categories_)
    :lambdas(lambdas_), alpha(alpha_), beta(beta_), probes(probes_), cat_indices(cat_indices_),
    noOfCategories(lambdas[0].size()),categories(categories_) {}

  double operator()(varray & kappa, double kappa2, int j/*, int cat*/)
  {
    double lam = 0.0;
    for(int i=0; i < noOfCategories; i++) {
      lam += cat_indices[i].size() * lambdas[j][i];
    }
    lam = lam/probes;


    double sum=0.0;
    for(int i=0;i < probes; ++i) {
      sum += (lam - lambdas[j][categories[i]]) * (lam - lambdas[j][categories[i]]);
    }

    return SAFE_EXP((kappa[j]-kappa2)*(.5*sum+beta)
        + (alpha+.5*probes-1)*log(kappa2/kappa[j]));
  }
};

class Phi_T{
private:
  double &phi;
  varray2d& MM;
  varray2d& S;
  varray2d& H;
  varray& tau;
  varray& beta;
  Random* rand;
  int probes;
  int samples;
  int *categories;
  int numberCategories;
public:
  Phi_T(double &phi_, varray2d& MM_, varray2d& S_, varray2d& H_, varray& tau_, varray& beta_,
	Random* rand_, int *categories_, int numberCategories_)
    : phi(phi_), MM(MM_), S(S_), H(H_), tau(tau_), beta(beta_), rand(rand_), probes(MM[0].size()), samples(MM.size()),
      categories(categories_), numberCategories(numberCategories_) {}

   void Update() {
    double sumMMmHTot=0;
    double sumSTot=0;
    for(int j=0; j<samples; ++j) {
      double sumMMmH=0;
      double sumS=0;
      for(int i=0; i<probes; ++i) {
        sumMMmH += S[j][i]*(MM[j][i]-H[j][i]-beta[j]);
        sumS += S[j][i]*S[j][i];
      }
      sumMMmHTot += sumMMmH*tau[j];
      sumSTot += sumS*tau[j];
    }
    phi=rand-> TruncNormal( sumMMmHTot/sumSTot, sqrt(1.0/(sumSTot)) );
  }

};

class Tau_T{
private:
  varray& tau;
  varray2d& PM;
  varray2d& MM;
  varray2d& S;
  varray2d& H;
  double &phi;
  varray& beta;
  double alpha;
  double p_beta;
  Random* rand;
  int probes;
  int samples;
public:
  Tau_T(varray& tau_, varray2d& PM_, varray2d& MM_, varray2d& S_, varray2d& H_, 
	double &phi_, varray& beta_, double alpha_, double p_beta_, Random* rand_)
    : tau(tau_), PM(PM_), MM(MM_), S(S_), H(H_), phi(phi_), beta(beta_), alpha(alpha_), 
      p_beta(p_beta_), rand(rand_), probes(PM[0].size()), samples(PM.size()) {}

  void Update()
  {
    for(int j=0; j<samples; ++j)
      {
	double sum=0;
	for(int i=0; i<probes; ++i)
	  {
	    sum += ( ( PM[j][i]-(S[j][i]+H[j][i]+beta[j]) ) * ( PM[j][i]-(S[j][i]+H[j][i]+beta[j]) ) + 
		 ( MM[j][i]-(phi*S[j][i]+H[j][i]+beta[j]) ) * ( MM[j][i]-(phi*S[j][i]+H[j][i]+beta[j]) ) );
	  }
	tau[j]=rand->Gamma( alpha+probes, p_beta+sum*.5 );
      }
  }
};

/* not in use - fixed by Emp Bayes like procedure
class A_T{
private:
  varray& a;
  varray2d& taug;
  varray& b;
  double p_mu;
  double p_tau;
  Random* rand;
  int genes;
  int conditions;
public:
  A_T(varray& a_, varray2d& taug_, varray& b_,
      double p_mu_, double p_tau_, Random* rand_)
    : a(a_), taug(taug_), b(b_),
      p_mu(p_mu_), p_tau(p_tau_), rand(rand_), genes(taug[0].size()), conditions(taug.size()) {}

  void Update()
  {
    for(int j=0; j<conditions; ++j)
      {
	double sum=0;
	for(int i=0; i<genes; ++i)
	  {
	    sum += log(taug[j][i]);
	  }
	a[j]=rand->Normal( (p_tau*p_mu+b[j]*sum)/(p_tau+b[j]*genes), 
			   sqrt( 1.0/(p_tau+b[j]*genes) ) );
      }
  }
};

class B_T{
private:
  varray& b;
  varray2d& taug;
  varray& a;
  double alpha;
  double beta;
  Random* rand;
  int genes;
  int conditions;
public:
  B_T(varray& b_, varray2d& taug_, varray& a_,
      double alpha_, double beta_, Random* rand_)
    : b(b_), taug(taug_), a(a_),
      alpha(alpha_), beta(beta_), rand(rand_), genes(taug[0].size()), conditions(taug.size()) {}

  void Update()
  {
    for(int j=0; j<conditions; ++j)
      {
	double rss=0;
	for(int i=0; i<genes; ++i)
	  {
	    rss += (log(taug[j][i])-a[j])*(log(taug[j][i])-a[j]);
	  }
	b[j]=rand->Gamma( alpha+genes*.5, beta+rss*.5 );
      }
  }
};

class Beta_T{
private:
  varray& beta;
  varray2d& PM;
  varray2d& MM;
  varray2d& S;
  varray2d& H;
  double& phi;
  varray& tau;
  double mu;
  double sigma;
  Random* rand;
  int probes;
  int samples;
public:
  Beta_T(varray& beta_, varray2d& PM_, varray2d& MM_, varray2d& S_, varray2d& H_,
      double& phi_, varray& tau_, double mu_, double sigma_, Random* rand_)
    : beta(beta_), PM(PM_),  MM(MM_), S(S_), H(H_), phi(phi_), tau(tau_),
      mu(mu_), sigma(sigma_), rand(rand_), probes(PM[0].size()), samples(PM.size()) {}

  void Update()
  {
    for(int j=1; j<samples; ++j)
      {
	double sum=0;
	for(int i=0; i<probes; ++i)
	  {
	    sum += PM[j][i]-S[j][i]-2*H[j][i]+MM[j][i]-phi*S[j][i];
	  }
	beta[j]=rand->Normal( tau[j]*sum/(probes*tau[j]+1/(2*sigma)), sqrt(1.0/(probes*tau[j]+1/(2*sigma))) );
      }
  }
};
*/

class MissingProbeSeqs {
  private:
  int *categories;
  int numberOfCategories;
  int probes,samples;
  int *unknownSeqProbePos;
  int numberOfUnknownSeqs;
  double *catDistribution;
  int *distinctCategories;
  double *gammas;
  double *discreteSamplingDistribution;
  double *cumulativeDiscreteSamplingDistribution;
  varray2d& H; 
  varray2d& lambdas;
  varray& eta;
  Random* rand;

  public:
  MissingProbeSeqs(int *categories_, int numberOfCategories_,int probes_, int samples_,
    int *unknownSeqProbePos_, int numberOfUnknownSeqs_, varray2d& H_, varray2d& lambdas_, varray& eta_,
    Random* rand_) :
    categories(categories_), numberOfCategories(numberOfCategories_), probes(probes_),samples(samples_),
    unknownSeqProbePos(unknownSeqProbePos_),numberOfUnknownSeqs(numberOfUnknownSeqs_),
    H(H_), lambdas(lambdas_), eta(eta_), rand(rand_) {

    catDistribution = new double[numberOfCategories];
    for(int i=0;i<numberOfCategories;i++)catDistribution[i]=0.0;
    int unk=0;
    for(int i=0; i <probes; i++) {
      if(numberOfUnknownSeqs && i==unknownSeqProbePos[unk] && unk < (numberOfUnknownSeqs-1)) unk++;
      else catDistribution[categories[i]] += 1.0;
    }

    for(int i=0;i<numberOfCategories;i++) catDistribution[i] /= (probes-numberOfUnknownSeqs);
    distinctCategories = new int[numberOfCategories];
    for(int i=0; i < numberOfCategories; i++) distinctCategories[i] = i;
    gammas = new double[numberOfCategories];
    discreteSamplingDistribution = new double[numberOfCategories];
    cumulativeDiscreteSamplingDistribution = new double[numberOfCategories];
  }

  ~MissingProbeSeqs(){
    delete[] catDistribution;
    delete[] distinctCategories;
    delete[] gammas;
    delete[] discreteSamplingDistribution; 
    delete[] cumulativeDiscreteSamplingDistribution;
  }

  void Update(){
    for(int i=0;i < numberOfUnknownSeqs; i++) {
      double ming=numeric_limits<double>::max();
      for(int k=0; k < numberOfCategories; k++) {
        gammas[k]=0;
        for(int j=0; j < samples; j++) {
          gammas[k] += 0.5*log(eta[j]) +
            log(pnorm(sqrt(eta[j])*lambdas[j][k])) +
            eta[j]*.5*(log(H[j][unknownSeqProbePos[i]]+1)-lambdas[j][k])*
              (log(H[j][unknownSeqProbePos[i]]+1)-lambdas[j][k]);
        }
        if(gammas[k] < ming) ming = gammas[k];
      }
      double sumOfDiscreteSamplingDist=0.0;
      for(int k=0; k < numberOfCategories; k++) {
        discreteSamplingDistribution[k] = catDistribution[k]*SAFE_EXP(-gammas[k]+ming);
        sumOfDiscreteSamplingDist += discreteSamplingDistribution[k];
      }
      for(int k=0; k < numberOfCategories; k++) {
        discreteSamplingDistribution[k] /= sumOfDiscreteSamplingDist;
      }

      cumulativeDiscreteSamplingDistribution[0] = discreteSamplingDistribution[0];
      for(int k=1; k < numberOfCategories; k++) {
        cumulativeDiscreteSamplingDistribution[k] = cumulativeDiscreteSamplingDistribution[k-1] + discreteSamplingDistribution[k];
      }

      categories[unknownSeqProbePos[i]] =
        rand->Discrete(cumulativeDiscreteSamplingDistribution,distinctCategories,numberOfCategories);

    }
  }
};

#endif // _BGX_UPDATES_HH
