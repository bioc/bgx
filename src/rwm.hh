/*
 *  This file is part of BGX, the Bayesian Gene eXpression program.
 *  Copyright (c) 2003-2004  Graeme Ambler <graeme@ambler.me.uk>
 *  Copyright (c) 2006 Ernest Turro <ernest.turro@ic.ac.uk>
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

#ifndef LOCAL_RWM_INCLUDED
#define LOCAL_RWM_INCLUDED

#include "rand.hh"

/* 1-dimensional Param */
template<class Accept, typename Update_T=std::valarray<double>, class Generator=Random>
class RWM{
 private:
  Update_T& Param;
  Accept& AccProb;
  Update_T& jmp_sizes;
  int batch_size;
  double optimal_AR;
  double start_jmp;

  Generator* rand;
  Update_T nAccept;
  int batch;
  Update_T nAcceptThisBatch;

  int length;
 public:
  RWM(Update_T& Param_, Accept& AccProb_, Update_T& jmp_sizes_, int batch_size_, double optimal_AR_, double start_jmp_, Generator* rand_)
    : Param(Param_), AccProb(AccProb_), jmp_sizes(jmp_sizes_), batch_size(batch_size_), optimal_AR(optimal_AR_),start_jmp(start_jmp_), rand(rand_),nAccept(Param),batch(1),nAcceptThisBatch(Param),length(Param.size())
    
  { nAccept=nAcceptThisBatch=0; 
  }

  void Update()
  {
    for(int i=0; i<length; ++i)
      {
	double proposal=rand->Normal(Param[i],jmp_sizes[i]);
	if(rand->Uniform() < AccProb(Param,proposal,i))
	  {
	    Param[i]=proposal;
	    ++nAccept[i];
	    ++nAcceptThisBatch[i];
	  }
      }
  }

  void Update_jmps() {
    double delta = std::min(0.01,pow(batch,-0.5));
    for(int i=0; i < length; i++) {
      if((double(nAcceptThisBatch[i])/batch_size) < optimal_AR) {
        jmp_sizes[i] *= exp(-delta);
      } else if((double(nAcceptThisBatch[i])/batch_size) > optimal_AR) {
        jmp_sizes[i] *= exp(delta);
      }
      nAcceptThisBatch[i]=0;
    }
    batch++;
  }

  void Reset() { for(int i=0; i<length; ++i) {nAccept[i]=nAcceptThisBatch[i]=0; }}

  void pAccept(int nTry, double* pAcc)
  {
    for(int i=0; i<length; ++i)
      pAcc[i] = nAccept[i]/(double) nTry;
  }
};


/* 2-dimensional Param */
template<class Accept, class Generator>
class RWM<Accept,std::valarray<std::valarray<double> >,Generator>{
 private:
  std::valarray<std::valarray<double> >& Param;
  Accept& AccProb;
//  double jmp_size;
  std::valarray<std::valarray<double> >& jmp_sizes;
  int batch_size;
  double optimal_AR;
  double start_jmp;
  Generator* rand;
  std::valarray<std::valarray<double> > nAccept;
  int batch;
  std::valarray<std::valarray<double> > nAcceptThisBatch;

  int length;
 public:
  RWM(std::valarray<std::valarray<double> >& Param_, Accept& AccProb_,
      std::valarray<std::valarray<double> >& jmp_sizes_, int batch_size_, double optimal_AR_, double start_jmp_, Generator* rand_)
    : Param(Param_), AccProb(AccProb_), jmp_sizes(jmp_sizes_), batch_size(batch_size_), optimal_AR(optimal_AR_),start_jmp(start_jmp_), rand(rand_),nAccept(Param),batch(1),nAcceptThisBatch(Param), length(Param.size())
  {
    for(int i=0; i<length; ++i) {
      nAccept[i].resize(Param[i].size());
      nAcceptThisBatch[i].resize(Param[i].size());
      nAccept[i]=0;
      nAcceptThisBatch[i]=0;
    }
  }

  void Update() {
  for(int j=0; j<length; ++j) {
	    for(unsigned int i=0; i<Param[j].size(); ++i) {
	      double proposal=rand->Normal(Param[j][i],jmp_sizes[j][i]);
	      if(rand->Uniform() < AccProb(Param,proposal,j,i)) {
          Param[j][i]=proposal;
          ++nAccept[j][i];
          ++nAcceptThisBatch[j][i];
        }
	    }
    }
  }

  void Update_jmps() {
    double delta = std::min(0.01,pow(batch,-0.5));
    for(int j=0; j<length; j++) {
      for(unsigned int i=0; i < Param[j].size(); i++) {
        if((double(nAcceptThisBatch[j][i])/batch_size) < optimal_AR) {
          jmp_sizes[j][i] *= exp(-delta);
        } else if((double(nAcceptThisBatch[j][i])/batch_size) > optimal_AR) {
          jmp_sizes[j][i] *= exp(delta);
        }
        nAcceptThisBatch[j][i]=0;
      }
    }
    batch++;
  }

  void Reset()
  {
    for(int j=0; j<length; ++j)
      for(unsigned int i=0; i<Param[j].size(); ++i){
	nAccept[j][i]=0;
  nAcceptThisBatch[j][i]=0;
      }
  }

  void pAccept(int nTry, double* pAcc)
  {
    int k=0;
    for(int i=0; i<length; ++i)
      for(unsigned int j=0; j<Param[i].size(); ++j)
	pAcc[k++] = nAccept[i][j]/(double) nTry;
  }
};

/* global Param */
template<class Accept, class Generator>
class RWM<Accept,double,Generator>{
private:
  double& Param;
  Accept& AccProb;
  double jmp_size;
  Generator* rand;
  double nAccept;
public:
  RWM(double& Param_, Accept& AccProb_, double jmp_size_, Generator* rand_)
    : Param(Param_), AccProb(AccProb_), jmp_size(jmp_size_), rand(rand_), nAccept(Param)
  { nAccept=0; }

  void Update()
  {
    double proposal=rand->Normal(Param,jmp_size);
    if(rand->Uniform() < AccProb(Param,proposal))
      {
	Param=proposal;
	++nAccept;
      }
  }

  void Reset(){ nAccept=0; }

  void pAccept(int nTry, double* pAcc)
  {
    *pAcc = nAccept/(double) nTry;
  }
};

#endif // LOCAL_RWM_INCLUDED
