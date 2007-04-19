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

#ifndef LOCAL_QNORM_INCLUDED
#define LOCAL_QNORM_INCLUDED

#ifdef __cplusplus
extern "C"{
#endif

// A really accurate qnorm routine using the "PPND16" algorithm of:
//  Wichura, M.J. (1988).
//   The Percentage Points of the Normal Distribution (Algorithm AS 241).
//   Applied Statistics 37, pp. 477-484.
double qnorm(double p);

// A slow and dirty qnorm routine I found on the web somewhere, accurate 
// to within 1 part in 2000 in the range 0.00001 to 0.99999 (and probably 
// outside that range too, but I've only tested it on [0.00001,0.99999]).
// Its only merits are that it can be written in very few lines of code 
// and that it only requires a single temporary variable.
double qnorm2(double p);

#ifdef __cplusplus
}
#endif

#endif // LOCAL_QNORM_INCLUDED
