/*
 *  This file is part of BGX, the Bayesian Gene eXpression program.
 *  Copyright (c) 2003-2004  Graeme Ambler <graeme@ambler.me.uk>
 *                2004       Anne-Mette Hein <a.hein@imperial.ac.uk>
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

#define BGXDEBUG 0

#ifdef USING_R
  extern "C" {
  #include <R.h> // for flushing console
  #include <R_ext/Utils.h> // we need to allow user interrupts from R
  #include <R_ext/Print.h> // R for windows doesn't like cout or printf, need Rprintf
  #if ( defined(HAVE_AQUA) || defined(Win32) )
    #define FLUSHCONSOLE {R_FlushConsole(); R_ProcessEvents();}
  #else
    #define FLUSHCONSOLE R_FlushConsole();
  #endif
  }
#endif

#include "rand.hh"
#include "rundir.hh"
#include "qnorm.h"
#include "sokal.hh"
#include "bgx_updates.hh"
#include <string>
#include <iostream>
#include <fstream>
#include <valarray>
#include <vector>
#include <time.h>
#include "bgx.hh"

#ifdef USING_R  /* Register routines, allocate resources. */
  #define R_NO_REMAP
  extern "C"{
  #include <Rinternals.h>
  #include <R_ext/Rdynload.h>
  }
  static const R_CMethodDef cMethods[] = {
    {"bgx", (DL_FUNC)&bgx, 34},
    {"freeBGXMemory", (DL_FUNC)&freeBGXMemory, 2},
    {NULL, NULL, 0}
  };
  void R_init_bgx(DllInfo *info) {
    R_registerRoutines(info, cMethods,NULL,NULL,NULL);
  }
  void R_unload_bgx(DllInfo *info) { } 
#endif

using namespace std;

enum { MINIMAL, BGXTRACE, ALL};

typedef valarray<double> array;
typedef valarray<valarray<double> > array2d;

void stringcpy(char* p, const string& s) {
  s.copy(p,string::npos);
  p[s.length()] = 0; // Add C-style string terminator
}

/* inline double TN_quantile(double p, double mu, double sigma) {
  return mu-sigma*qnorm(pnorm(mu/sigma)*(1-p));
} */


// The following are global so that we can delete them in a separate function if there is a user interrupt in R
static double *xave, *yave, *muave, *sacc, *hacc, *muacc, *sigmaacc, *lambdaacc, *etaacc;
static ofstream *sigma_, *lambda_, *XS, *YS;
static fstream *mu_;
  // Metropolis/Gibbs objects
static  S_T *AccS;
static  RWM<S_T,array2d> *S;
static  H_T *AccH;
static  RWM<H_T,array2d> *H;
static  Mu_T *AccMu;
static  RWM<Mu_T,array2d> *Mu;
static  Sigma_T *AccSigma;
static  RWM<Sigma_T,array2d> *Sigma;
static  Lambda_T *AccLambda;
static  RWM<Lambda_T,array2d> *Lambda;
static  Eta_T *AccEta;
static  RWM<Eta_T> *Eta;
static  Phi_T *Phi;
static  Tau_T *Tau;

//extern "C"
//void freeBGXMemory(int *output, int *numberGenesToWatch);

void bgx(double* pm, double* mm, int* samples, int* conditions, 
	 int* probes, int* genes, int *numberCategories, int* numberGenesToWatch, int* samplesets, int* probesets,
   int* categories, int* unknownProbeSeqs, int* numberOfUnknownProbeSeqs, int* whichGenesToWatch, 
   int* whichProbesToWatch, int* iter, int* burnin, double* s_jmp, double* h_jmp,
	 double* mu_jmp, double* sigma_jmp, double* lambda_jmp, 
	 double* eta_jmp, bool* adaptive, int* batch_size, double* optimalAR, 
	 int* subsample, int* output, char** dirname, char **basepath, int* seed, char **sampleNames, string * affinityPlotFile=NULL, string * geneNamesFile=NULL) {

  // SETUP. declaration/initialisation of some variables up here for common code structure with parallel BGX
  string run_dir;

  array2d s(*samples),h(*samples),PM(*samples),MM(*samples), lambda(*samples),
    mu(*conditions),sigma(*conditions), s_jmps(*samples), h_jmps(*samples), 
    mu_jmps(*conditions), sigma_jmps(*conditions), lambda_jmps(*samples);

  for(int j=0; j < *samples; ++j) {
    s[j].resize(*probes);
    s_jmps[j].resize(*probes);
    h[j].resize(*probes);
    h_jmps[j].resize(*probes);
    MM[j].resize(*probes);
    PM[j].resize(*probes);
    lambda[j].resize(*numberCategories);
    lambda_jmps[j].resize(*numberCategories);
    for(int i=0; i < *probes; i++) {
      s_jmps[j][i]=*s_jmp;
      h_jmps[j][i]=*h_jmp;
    }
    for(int k=0; k < *numberCategories; k++) lambda_jmps[j][k] = *lambda_jmp;
  }
  for(int j=0; j < *conditions; ++j) {
    mu[j].resize(*genes);
    mu_jmps[j].resize(*genes);
    sigma[j].resize(*genes);
    sigma_jmps[j].resize(*genes);
    for(int i=0; i<*genes; ++i){
      mu_jmps[j][i]=*mu_jmp;
      sigma_jmps[j][i]=*sigma_jmp;
    }
  }

  // Throughout the code we use eta, b, sigma  to refer to eta^{-2}, b^{-2}, sigma^{-2}
  array tau(*samples), beta(*samples), eta(*samples), eta_jmps(*samples), a(*conditions), b(*conditions);
  for(int j=0; j < *samples; j++) eta_jmps[j] = *eta_jmp;

  double phi;
  
  ofstream  muave_, xave_, yave_, sacc_, hacc_, muacc_, sigmaacc_, lambdaacc_, etaacc_,
    tau_, phi_,  eta_, unkcats_;

  // Metropolis/Gibbs objects
/*
  S_T *AccS;
  RWM<S_T,array2d> *S;
  H_T *AccH;
  RWM<H_T,array2d> *H;
  Mu_T *AccMu;
  RWM<Mu_T,array2d> *Mu;
  Sigma_T *AccSigma;
  RWM<Sigma_T,array2d> *Sigma;
  Lambda_T *AccLambda;
  RWM<Lambda_T,array2d> *Lambda;
  Eta_T *AccEta;
  RWM<Eta_T> *Eta;
  Phi_T *Phi;
  Tau_T *Tau;
*/

  // Initialise "averages"
  muave = new double[*genes**conditions];
  for(int i=0; i<*genes**conditions; ++i) muave[i]=0;

  if(*output >= BGXTRACE){
    xave = new double[*probes**samples];
    yave = new double[*probes**samples];
    for(int i=0; i<*probes**samples; ++i) xave[i]=yave[i]=0;
  }

  // END SETUP

  Random rand(*seed);

  // Create output file and directory
  char tmpstr[80];
  strcpy(tmpstr, *basepath); // Prepend basepath
  strcat(tmpstr,"/run");
  run_dir = rundir(tmpstr);

  // Copy affinity plot file to run_dir if available
  if(affinityPlotFile!=NULL) {
    ifstream in(affinityPlotFile->c_str());
    unsigned long int lastslash = affinityPlotFile->find_last_of('/');
    if(lastslash != string::npos) *affinityPlotFile = affinityPlotFile->substr(lastslash+1);;
    ofstream out((run_dir + '/' + *affinityPlotFile).c_str());
    out << in.rdbuf();
  }

 // Copy geneNames file to run_dir if available
  if(geneNamesFile!=NULL) {
    ifstream in(geneNamesFile->c_str());
    unsigned long int lastslash = geneNamesFile->find_last_of('/');
    if(lastslash != string::npos) *geneNamesFile = geneNamesFile->substr(lastslash+1);;
    ofstream out((run_dir + '/' + *geneNamesFile).c_str());
    out << in.rdbuf();
  }

  // Write out a summary of the call
  string filename = run_dir+"/summary.txt";
  ofstream summary(filename.c_str());
  summary << "The program was called with the following parameter values:" << endl << endl
	  << "Number of samples:\t" << *samples << endl
	  << "Number of conditions:\t" << *conditions << endl
	  << "Number of probes:\t" << *probes << endl
	  << "Number of genes:\t" << *genes << endl 
    << "Number of categories:\t" << *numberCategories << endl
    << "Number of probes with no sequence information:\t" << *numberOfUnknownProbeSeqs << endl
    << "Number of genes to monitor fully:\t" << *numberGenesToWatch << endl
	  << "Subsampling interval:\t" << *subsample << endl
	  << "Number of burn-in sweeps:\t" << *burnin << endl
	  << "Number of post burn-in sweeps:\t" << *iter << endl
	  << "Spread of jumps in S:\t" << *s_jmp << endl
	  << "Spread of jumps in H:\t" << *h_jmp << endl
	  << "Spread of jumps in mu:\t" << *mu_jmp << endl
	  << "Spread of jumps in sigma:\t" << *sigma_jmp << endl
	  << "Spread of jumps in lambda:\t" << *lambda_jmp << endl
	  << "Spread of jumps in eta:\t" << *eta_jmp << endl
	  << "Adaptive MCMC:\t" << *adaptive << endl
	  << "Batch size:\t" << *batch_size << endl
	  << "Optimal acceptance ratio:\t" << *optimalAR << endl
	  << "Name of output directory:\t" << run_dir << endl
	  << "Seed for random number generator:\t" << *seed << endl
    << "CEL files:\t";
    for(int s =0; s < *samples-1; s++) summary << sampleNames[s] << ", ";
    summary << sampleNames[*samples-1] << endl;

  // In R, dirname is set to the ouput directory
  stringcpy(*dirname,run_dir);
  
  // Open outoput files
  filename=run_dir+"/muave";
  muave_.open(filename.c_str());
  mu_ = new fstream[*conditions];
  string num;
  for(int j=0; j<*conditions; ++j){
    int_to_string(j+1,num);
    filename=run_dir+"/mu."+num;
    // Open for in+out because we want to read these again later to calculate MCSE
    mu_[j].open(filename.c_str(), fstream::in | fstream::out | fstream::trunc);
  }
    
  if(*output >= BGXTRACE) {
    filename=run_dir+"/xave";
    xave_.open(filename.c_str());
    filename=run_dir+"/yave";
    yave_.open(filename.c_str());
    sigma_ = new ofstream[*conditions];
    for(int j=0; j<*conditions; ++j){
      int_to_string(j+1,num);
      filename=run_dir+"/sigma2."+num;
      sigma_[j].open(filename.c_str());
    }

    lambda_ = new ofstream[*samples];
    filename=run_dir+"/tau2";
    tau_.open(filename.c_str());
    filename=run_dir+"/phi";
    phi_.open(filename.c_str());

    for(int i=0; i < *samples; i++) {
      int_to_string(i+1,num);
      filename=run_dir+"/lambda."+num;
      lambda_[i].open(filename.c_str());
    }
    filename=run_dir+"/eta2";
    eta_.open(filename.c_str());


    if(*numberOfUnknownProbeSeqs) {
      filename=run_dir+"/unknownCategories";
      unkcats_.open(filename.c_str());
    }
  }

  if(*output >= ALL) {
    // Acceptance probabilities
    sacc = new double[*probes**samples];
    hacc = new double[*probes**samples];
    muacc = new double[*genes**conditions];
    sigmaacc = new double[*genes**conditions];
    etaacc = new double[*samples]; for(int i=0; i < *samples; i++) etaacc[i] = 0;
    lambdaacc = new double[*samples**numberCategories];
    for(int i=0; i<*probes**samples; ++i)
      sacc[i]=hacc[i]=0;
    for(int i=0; i<*genes**conditions; ++i)
      muacc[i]=sigmaacc[i]=0;
    for(int i=0; i<*samples**numberCategories; ++i)
      lambdaacc[i]=0;

    // Open output files
    filename=run_dir+"/sacc";
    sacc_.open(filename.c_str());
    filename=run_dir+"/hacc";
    hacc_.open(filename.c_str());
    filename=run_dir+"/muacc";
    muacc_.open(filename.c_str());
    filename=run_dir+"/sigmaacc";
    sigmaacc_.open(filename.c_str());
    filename=run_dir+"/lambdaacc";
    lambdaacc_.open(filename.c_str());
    filename=run_dir+"/etaacc";
    etaacc_.open(filename.c_str());
  }

  if(*numberGenesToWatch>0){
    XS = new ofstream[*samples];
    YS = new ofstream[*samples];
    for(int j=0; j<*samples; ++j){
      int_to_string(j+1,num);
      filename=run_dir+"/x."+num;
      XS[j].open(filename.c_str());
      filename=run_dir+"/y."+num;
      YS[j].open(filename.c_str());
    }
  }

  // Read data into C++-style arrays for ease of manipulation
  for(int j=0; j<*samples; ++j){
    for(int i=0; i<*probes; ++i){
      PM[j][i]=*(pm+i+*probes*j);
      MM[j][i]=*(mm+i+*probes*j);
    }
  }

  // Define parameters and pick some vaguely sensible values for 
  // the initial values.

  phi = 0.18;
  for(int j=0; j<*samples; ++j){
    tau[j]=0.0005/**/*rand.Uniform(0.5,2)/**/;
    beta[j]=0;
    double mean_lh=0; double var_lh=0;
    const double tmp2=log(1.1/**/+rand.Uniform(0,1)/**/);
    for(int i=0; i<*probes; ++i){
      double tmp = (PM[j][i]-MM[j][i])/(1-phi);
      s[j][i] = 0.1>tmp ? 0.1 : tmp; 
      tmp = (phi*PM[j][i]-MM[j][i])/(phi-1);
      if(0.1>tmp){
        h[j][i]=0.1;
        mean_lh += tmp2;
      }else{
        h[j][i]=tmp;
        mean_lh+=log(tmp+1);
      }
    }
    mean_lh /= *probes;
    for(int i=0; i<*probes; ++i){
      double L = log(h[j][i]+1);
      var_lh += (L-mean_lh)*(L-mean_lh);
    }
    var_lh /= *probes-1;
    for(int i=0; i < *numberCategories; i++) {
      lambda[j][i]=mean_lh;
      //etas[j][i]=1.0/var_lh;
    }
    eta[j]=1.0/var_lh;
  }
  
  // AMs addition: Obtaining Empirical Bayes like
  // estiamtes of a and b^2:

  array2d *ebS = new array2d(*samples);
  // setting empirical Bayes values of signals: (*ebS)
  vector<pair<int,int> > badSets; // probe sets which have PM<=MM for all probes

  for(int j=0; j<*samples; ++j){
    (*ebS)[j].resize(*probes);  
    int probeCounter = 0;
    for(int i=0; i<*genes; ++i){        
      // If PM>MM: (*ebS) = PM-MM. Furthermore, for each probe set: 
      // find minimum positive diff and number probe pairs with PM>MM. 
      int positive = 0;
      // We test minimumDiff == HUGE_VAL below so be careful changing this!
      double minimumDiff = HUGE_VAL;
      for(int p=0; p<probesets[i]; ++p){
        if ( PM[j][probeCounter+p] > MM[j][probeCounter+p] ){
          (*ebS)[j][probeCounter+p] = PM[j][probeCounter+p] - MM[j][probeCounter+p];
          positive++;
          if( minimumDiff > (*ebS)[j][probeCounter+p] ){
            minimumDiff = (*ebS)[j][probeCounter+p];
          }
        }
      }
      if(minimumDiff == HUGE_VAL){ // PM<=MM for all probes in probeset
        pair<int,int> tempPair(j,i);
        badSets.push_back(tempPair);
      }
      for(int p=0; p<probesets[i]; ++p){
        if(minimumDiff==HUGE_VAL || PM[j][probeCounter+p]==MM[j][probeCounter+p]) (*ebS)[j][probeCounter+p] = 0.0;
        // If PM<MM: (*ebS) = fraction of minimum positive PM-MM
        else if ( PM[j][probeCounter+p] < MM[j][probeCounter+p] ){    
          (*ebS)[j][probeCounter+p] = (double(positive)/double(probesets[i]))*minimumDiff;
        }
      }
      probeCounter += probesets[i];
    }
  }

  // For each (g,c), calculate variance over set consisting of
  // the values log((*ebS)_{g,j,c,r}+1), j=1,...,J_g, r=1,...,R_c

  bool badData = false;
  {
    array2d vars(*conditions);
    int sampleCounter = 0;
    for(int j=0; j<*conditions; ++j){
      vars[j].resize(*genes);
      int probeCounter = 0;
      for(int i=0; i<*genes; ++i){
        double tempS = 0;
        double tempSS = 0;                
          
      for(int sam=0; sam<samplesets[j]; ++sam){
        for(int p=0; p<probesets[i]; ++p){
          double logEbS=log((*ebS)[sampleCounter+sam][probeCounter+p]+1.0);
            tempS += logEbS;
            tempSS += logEbS*logEbS;
          }
        }
        probeCounter += probesets[i];
        vars[j][i] = (1.0/(double(probesets[i] * samplesets[j])-1.0))*
                     (tempSS-(tempS*tempS/double(probesets[i]*samplesets[j])));
       }
       sampleCounter += samplesets[j];

       // For each condition, calculate the mean and the variance
       // (or precision) of the log-variances log(vars[j][i]):
       double tempS = 0.0;
       double tempSS = 0.0;
       double minVar = HUGE_VAL;
       for(int i=0; i<*genes; ++i)
         if(vars[j][i]<minVar && vars[j][i]>0) minVar=vars[j][i];
       if(minVar == HUGE_VAL){ // PM<=MM for all probe sets.
         // Set a,b to default values.
         badData = true;
         a[j] = 0;
         b[j] = 3;
       } else {
         minVar=log(minVar);
         for(int i=0; i<*genes; ++i){
           double logVar = 0.0;
           if(vars[j][i]>0.0) logVar = log(vars[j][i]);
           else logVar = minVar;
           tempS += logVar;
           tempSS += logVar*logVar;
        }
        a[j] = tempS / double(*genes);
        // NB! b[j] is the precision (=1/variance): b^{-2}
        b[j] = (double(*genes)-1)/(tempSS-(tempS*tempS / double(*genes)));
      }
      //cout << "empirical Bayes-like estimates of a[" << j << "]: " << a[j]
        //   << "  b[" << j << "]: " << b[j] << endl;
    }
  }
  delete ebS;

  filename=run_dir+"/empBayesEst";
  ofstream empBayesEst_(filename.c_str());
  empBayesEst_ << "Empirical Bayes like estimates of a and b^{-2} for each condition:" 
               << endl << endl;
  empBayesEst_ << '\t' << "condition"  <<  '\t' << "a"  <<  '\t' << "b^{-2}"  << endl;
  for(int j=0; j<*conditions; j++){
    empBayesEst_ << '\t'  << j << '\t'  <<  a[j] << '\t'  <<  b[j] << endl;
  }
  if(badSets.size()>1){
#ifndef USING_R
    printf("Warning: There were %d probe sets with PM<=MM for all probes.\n", (int)badSets.size());
#else
    Rprintf((char *)"Warning: There were %d probe sets with PM<=MM for all probes.\n", (int)badSets.size());
#endif
    empBayesEst_ << endl << "There were " << badSets.size() 
         << " probe sets with PM<=MM for all probes." << endl;
  } else if(badSets.size()==1){
#ifndef USING_R
    printf("Warning: There was %d probe set with PM<=MM for all probes.\n", (int)badSets.size());
#else
    Rprintf((char *)"Warning: There was %d probe set with PM<=MM for all probes.\n", (int)badSets.size());
#endif
    empBayesEst_ << endl << "There was " << badSets.size() 
                 << " probe set with PM<=MM for all probes." << endl;
  }
  empBayesEst_.close();
             
  // AMs addition end!

  // Setting initial values for mu_gc and sigma^2_gc
  if(badData){ // PM<=MM for all probe sets.
    for(int i=0; i<*genes; ++i){
      for(int j=0; j<*conditions; ++j){
        mu[j][i] = 0;
        sigma[j][i] = 1;
      }
    }
  } else {
    int p_loc=0;
    double minVar = HUGE_VAL;

    for(int i=0; i<*genes; ++i){
      int s_loc=0;

      for(int j=0; j<*conditions; ++j){
        double mean_ls=0; double var_ls=0;
        for(int sam=0; sam<samplesets[j]; ++sam){
	        for(int p=0; p<probesets[i]; ++p)
	          mean_ls += log(s[sam+s_loc][p+p_loc]+1);
        }
        mean_ls /= samplesets[j]*probesets[i];
        for(int sam=0; sam<samplesets[j]; ++sam){
          for(int p=0; p<probesets[i]; ++p){
            double L = log(s[sam+s_loc][p+p_loc]+1);
            var_ls += (L-mean_ls)*(L-mean_ls);
          }
        }
        var_ls /= samplesets[j]*probesets[i]-1;
        mu[j][i]=mean_ls;
        if(var_ls>0.0){
          sigma[j][i]=1/var_ls;
          if(var_ls<minVar) minVar=var_ls;
        }else{ // AMs addititon begin
          sigma[j][i]=1/minVar;
        }  // AMs addition end
        s_loc += samplesets[j];
      }
      p_loc += probesets[i];
    }
  }

  // Create a 2d vector with the category indices
  int unk=0;
  vector<vector<int> > cat_indices(*numberCategories);
  for(int i =0; i < *probes; i++) {
    if(*numberOfUnknownProbeSeqs && i==unknownProbeSeqs[unk] && unk < (*numberOfUnknownProbeSeqs-1)) unk++;
    else cat_indices[categories[i]].push_back(i);
  }

  // This Gibbs object is for updating categories for probes with unknown sequences
  MissingProbeSeqs I(categories,*numberCategories, *probes,*samples, 
    unknownProbeSeqs, *numberOfUnknownProbeSeqs, h, lambda, eta, &rand);

  // Initialise the update objects. For random walk Metropolis (RWM) steps, this 
  // involves creating an acceptance probability functor (e.g. S_T AccS(...)) 
  // and a RWM object (e.g. RWM<S_T,array2d> S(...)). For Gibbs moves we just
  // create a single update object (e.g. A_T A(...)).
  // These are all defined in bgx_updates.hh.

  AccS = new S_T(PM,MM,h,phi,mu,sigma,tau,beta,probesets,samplesets,categories );
  S = new RWM<S_T,array2d>(s,*AccS,s_jmps,*batch_size,*optimalAR,*s_jmp,&rand);

  AccH = new H_T(PM,MM,s,phi,lambda,eta,tau,beta,categories);
  H = new RWM<H_T,array2d>(h,*AccH,h_jmps,*batch_size,*optimalAR, *h_jmp, &rand);

  AccMu = new Mu_T(s,sigma,probesets,samplesets);
  Mu = new RWM<Mu_T,array2d>(mu,*AccMu,mu_jmps,*batch_size,*optimalAR,*mu_jmp,&rand);

  AccSigma = new Sigma_T(s,mu,a,b,probesets,samplesets);
  Sigma = new RWM<Sigma_T,array2d>(sigma,*AccSigma,sigma_jmps,*batch_size,*optimalAR,*sigma_jmp,&rand);

  AccLambda = new Lambda_T(h,eta,0,0.001, cat_indices);
  Lambda = new RWM<Lambda_T,array2d>(lambda,*AccLambda,lambda_jmps,*batch_size,*optimalAR,*lambda_jmp,&rand);

  AccEta = new Eta_T(h,lambda,0.001,0.001,cat_indices,categories);
  Eta = new RWM<Eta_T>(eta,*AccEta,eta_jmps,*batch_size,*optimalAR,*eta_jmp,&rand);

  Phi = new Phi_T(phi,MM,s,h,tau,beta,&rand, categories, *numberCategories);

  Tau = new Tau_T(tau,PM,MM,s,h,phi,beta,0.001,0.001,&rand);

#if BGXDEBUG
  ofstream fullMuTrace_((run_dir + "/fullMuTrace").c_str());
#endif

  time_t start_time, end_time;
  time(&start_time);

  // Run the Markov chain: Burn-in phase
  for(int sweep=0; sweep<*burnin; ++sweep){
#ifndef USING_R
    printf("Burn-in sweep %d of %d\r", sweep , *burnin);
    if(sweep%16 == 0) fflush(stdout);
#else
    Rprintf((char *)"Burn-in sweep %d of %d\r", sweep , *burnin);
    if(sweep%16 == 0) FLUSHCONSOLE
    R_CheckUserInterrupt();
#endif

#if BGXDEBUG
  if(sweep%8 == 0) {
    for(int c=0; c < *conditions; c++) {
      for(int g=0; g < *genes; g+=40) {
        fullMuTrace_ << mu[c][g] << " ";
      }
      fullMuTrace_ << endl;
    }
  }
#endif

    if(*numberOfUnknownProbeSeqs) I.Update();
      
    S->Update();
    H->Update();

    Mu->Update();
    Sigma->Update();

    Lambda->Update();
    Eta->Update();

    Tau->Update();
    Phi->Update();

    if(*adaptive && (sweep+1)%*batch_size==0){
      S->Update_jmps();
      H->Update_jmps();
      Mu->Update_jmps();
      Sigma->Update_jmps();
      Lambda->Update_jmps();

      Eta->Update_jmps();
    }
    
  }
#ifndef USING_R
  printf("%d burnin sweeps completed.\n", *burnin);
#else
  Rprintf((char *)"%d burnin sweeps completed.\n", *burnin);
#endif

  S->Reset();
  H->Reset();

  Mu->Reset();
  Sigma->Reset();

  Lambda->Reset();
  Eta->Reset();


  // Run the Markov chain: Sampling phase
  for(int sweep=0; sweep<*iter; ++sweep){
#ifndef USING_R
    printf("Post burn-in sweep %d of %d\r", sweep, *iter);
    if(sweep%16 == 0) fflush(stdout);
#else
    Rprintf((char *)"Post burn-in sweep %d of %d\r", sweep, *iter);
    if(sweep%16 == 0) FLUSHCONSOLE
    R_CheckUserInterrupt();
#endif

#if BGXDEBUG
  if(sweep%8 == 0) {
    for(int c=0; c < *conditions; c++) {
      for(int g=0; g < *genes; g++) {
        fullMuTrace_ << mu[c][g] << " ";
      }
      fullMuTrace_ << endl;
    }
  }
#endif

    if(*numberOfUnknownProbeSeqs) I.Update();

    S->Update();
    H->Update();

    Mu->Update();
    Sigma->Update();

    Lambda->Update();
    Eta->Update();

    Tau->Update();
    Phi->Update();


    // Not adapting jumps in post burn-in phase
/*    if(*adaptive && (sweep+1)%*batch_size==0){
      S->Update_jmps();
      H->Update_jmps();
      Mu->Update_jmps();
      Sigma->Update_jmps();
      Lambda->Update_jmps();

      Eta->Update_jmps();
    } */


    /* EACH SUBSAMPLE */
    if(!(sweep%*subsample)){
      for(int c=0; c<*conditions; ++c){
        // Update means of mus and sigmas
        for(int g=0; g<*genes; ++g){
          muave[g+*genes*c] += mu[c][g];
        }
      }

      
      for(int c=0; c<*conditions; ++c){
        for(int g=0; g<*genes; ++g) {
          mu_[c] << mu[c][g] << '\t';
        }
        mu_[c] << '\n';
      }

      if(*output >= BGXTRACE){
        for(int sample=0; sample<*samples; ++sample){
          // Update means of S and H
          for(int p=0; p<*probes; ++p){
            xave[p+*probes*sample] += log(s[sample][p]+1);
            yave[p+*probes*sample] += log(h[sample][p]+1);
          }
        }

        for(int c=0; c<*conditions; ++c){
          for(int g=0; g<*genes; ++g) {
            sigma_[c] << 1.0/sigma[c][g] << '\t';
          }
          sigma_[c] << '\n';
        }
        
        phi_ << phi << endl;

        for(int sample=0; sample<*samples; ++sample){
          tau_ << 1/tau[sample] << '\t';
          eta_ << 1/eta[sample] << '\t';
          for(int c=0; c<*numberCategories; ++c){
            lambda_[sample] << lambda[sample][c] << '\t';
           // etas_[sample] << 1.0/etas[sample][c] << '\t';
          }
          lambda_[sample] << '\n';
         // etas_[sample] << '\n';
        }
        tau_ << '\n';
     // lambda_ << '\n';
        eta_ << '\n';

        if(*numberOfUnknownProbeSeqs) {
          for(int i=0; i < *numberOfUnknownProbeSeqs;i++) unkcats_ << categories[unknownProbeSeqs[i]] << '\t';
          unkcats_ << endl;
        }
      }
 
      if(*numberGenesToWatch>0) {
        for(int sample=0; sample<*samples; ++sample){
          for(int g=0; g<*numberGenesToWatch; g++){
            for(int p=0; p<probesets[whichGenesToWatch[g]]; p++){
              // output x and y
              XS[sample] << log(s[sample][(whichProbesToWatch[g]-1)+p]+1) << '\t';
              YS[sample] << log(h[sample][(whichProbesToWatch[g]-1)+p]+1) << '\t';
            }
          }
          XS[sample] << '\n';
          YS[sample] << '\n';
        }
      }
    }
  }
  // End of simulation

  time(&end_time);
  int duration = (int)difftime(end_time,start_time);

#ifndef USING_R
  printf("%d post burn-in sweeps completed.\n", *iter);
  printf("MCMC duration: %dh %dm %ds\n", duration/3600, duration%3600/60, duration%60);
#else
  Rprintf((char *)"%d post burn-in sweeps completed.\n", *iter);
  Rprintf((char *)"MCMC duration: %dh %dm %ds\n", duration/3600, duration%3600/60, duration%60);
#endif

  summary << "MCMC duration:\t" <<  duration/3600 << "h " << duration%3600/60 << "m " << duration%60 << "s\n";
  summary.close();

#if BGXDEBUG
  fullMuTrace_.close();
#endif


  for(int c=0; c<*conditions; ++c){
    for(int g=0; g<*genes; ++g){
      muave_ << muave[g+*genes*c]**subsample/(*iter) << '\t';
    }
    muave_ << '\n';
  }
  muave_.close();

  if(*output >= BGXTRACE){
    for(int sample=0; sample<*samples; ++sample){
      for(int p=0; p<*probes; ++p){
        xave_ << xave[p+*probes*sample]**subsample/(*iter) << '\t';
        yave_ << yave[p+*probes*sample]**subsample/(*iter) << '\t';
      }
      xave_ << '\n';
      yave_ << '\n';
    }
    xave_.close();
    yave_.close();
  }


  if(*output >= BGXTRACE) {
    for(int c=0; c<*conditions; ++c) sigma_[c].close();
    tau_.close();
    phi_.close();
    for(int i=0; i < *samples; i++) lambda_[i].close();
    eta_.close();
    if(*numberOfUnknownProbeSeqs) unkcats_.close();
  }
  
  if(*output >= ALL) {
    // Calculate acceptance probabilities
    S->pAccept(*iter,sacc);
    H->pAccept(*iter,hacc);
    Mu->pAccept(*iter,muacc);
    Sigma->pAccept(*iter,sigmaacc);
    Lambda->pAccept(*iter,lambdaacc);
    Eta->pAccept(*iter,etaacc);
    for(int c=0; c<*conditions; ++c){
	    for(int g=0; g<*genes; ++g){
        muacc_ << muacc[g+*genes*c] << '\t';
        sigmaacc_ << sigmaacc[g+*genes*c] << '\t';
      }
      muacc_ << '\n';
      sigmaacc_ << '\n';
    }
  
    for(int sample=0; sample<*samples; ++sample){
	    for(int p=0; p<*probes; ++p){
        sacc_ << sacc[p+*probes*sample] << '\t';
        hacc_ << hacc[p+*probes*sample] << '\t';
      }
      sacc_ << '\n';
      hacc_ << '\n';
      for(int i=0; i < *numberCategories; i++) lambdaacc_ << lambdaacc[sample**numberCategories+i] << '\t';
      etaacc_ << etaacc[sample] << '\t';
      lambdaacc_ << '\n';
    }
    etaacc_ << '\n';
    
    sacc_.close();
    hacc_.close();
    muacc_.close();
    sigmaacc_.close();
    lambdaacc_.close();
    etaacc_.close();
  }

  if(*numberGenesToWatch>0) {
    for(int sample=0; sample<*samples; ++sample){
      XS[sample].close();
      YS[sample].close();
    }
  }

  // IACT and MCSE
  ofstream mcse_((run_dir + "/MCSE").c_str());
  ofstream iact_((run_dir + "/IACT").c_str());
  double * mcse = new double[*genes * *iter/(*subsample)];
  for(int c=0; c < *conditions; c++) {
    mu_[c].seekg(0, ios::beg);
    
    for(int i=0; i < *iter/(*subsample); i++) {
      for(int g=0; g < *genes; g++) mu_[c] >> mcse[g**iter/(*subsample) + i];
    }

    for(int g=0; g < *genes; g++) {
      double var=0;
      double tau=0;
      int m=0;
      if(sokal(*iter/(*subsample),&mcse[g**iter/(*subsample)],&var,&tau,&m)!=0) {
        mcse_ << *iter/(*subsample) << " ";
        iact_ << "NA" << " ";
      } else {
        mcse_ << tau*var**subsample/(*iter) << " ";
        iact_ << tau << " ";
      }
    }
    mcse_ << endl;
    iact_ << endl;
  }
  mcse_.close();
  iact_.close();
  for(int c=0; c<*conditions; ++c) mu_[c].close();
  delete [] mcse;

  // In the parallel version, memory that is allocated by all nodes is deallocated here rather than above.
  // We do the same here for code structure compatibility.

#ifndef USING_R
  freeBGXMemory(output, numberGenesToWatch);
#endif

#ifndef USING_R
  printf("Leaving bgx()...\n");
#else
  Rprintf((char *)"Leaving bgx()...\n");
#endif
}

void freeBGXMemory(int *output, int *numberGenesToWatch) {
#ifndef USING_R
  printf("Deallocating memory...\n");
#else
  Rprintf((char *)"Deallocating memory...\n");
#endif
  delete [] mu_;
  delete [] muave;
  delete AccS;
  delete S;
  delete AccH;
  delete H;
  delete AccMu;
  delete Mu;
  delete AccSigma;
  delete Sigma;
  delete AccLambda;
  delete Lambda;
  delete AccEta;
  delete Eta;
  delete Phi;
  delete Tau;

  if(*output >= BGXTRACE) {
    delete [] sigma_;
    delete [] lambda_;
    delete [] xave;
    delete [] yave;
  }

  if(*output >= ALL) {
    delete [] muacc;
    delete [] sigmaacc;
    delete [] sacc;
    delete [] hacc;
    delete [] lambdaacc;
    delete [] etaacc;
  }

  if(*numberGenesToWatch>0) {
    delete [] XS;
    delete [] YS;
  }
}
