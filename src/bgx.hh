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
#ifndef _BGX_HH
#define _BGX_HH
extern "C"
void bgx(double* pm, double* mm, int* samples, int* conditions, 
	  int* probes, int* genes, int *numberCategories, int* numberGenesToWatch, int* samplesets, int* probesets, 
    int* categories, int* unknownProbeSeqs, int* numberOfUnknownProbeSeqs, int* whichGenesToWatch, 
    int* whichProbesToWatch, int* iter, int* burnin, double* s_jmp, double* h_jmp, 
	  double* mu_jmp, double* tau_jmp, double* lambda_jmp, 
	  double* eta_jmp, bool* adaptive, int* batch_size, double* optimalAR,
	  int* subsample, int* output, char** dirname, char **basepath, int* seed, char ** sampleNames, std::string * plotInfoFile, std::string * geneNamesFile);


extern "C"
void freeBGXMemory(int *output, int *numberGenesToWatch);

#endif // _BGX_HH
