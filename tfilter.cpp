/*
 *  tfilter.cpp
 *  CrescendoPlayThru
 *
 *  Created by David McClain on 10/18/07.
 *  Copyright 2007 Refined Audiometrics Laboratory. All rights reserved.
 *
 */
/* -----------------------------------------------------------------------------
 Copyright (c) 2016 Refined Audiometrics Laboratory, LLC
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions
 are met:
 1. Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.
 2. Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.
 3. The names of the authors and contributors may not be used to endorse
 or promote products derived from this software without specific prior
 written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE AUTHORS AND CONTRIBUTORS ``AS IS'' AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 SUCH DAMAGE.
 ------------------------------------------------------------------------------- */

#include <math.h>
#include <memory.h>
#include "tfilter.h"
#include "useful_math.h"

// -------------------------------------------------------------------------------------
TFilter::TFilter(Float64 *pcoffs)
{
	memcpy(m_coffs, pcoffs, 5*sizeof(Float64));
	reset();
}

TFilter::~TFilter()
{}

void TFilter::reset()
{
	memset(m_state, 0, 2*sizeof(Float64));
}

#if 1
void TFilter::rdz()
{
    m_state[0] = ::rdz(m_state[0]);
    m_state[1] = ::rdz(m_state[1]);
}
#endif

// -------------------------------------------------------------------------------------
inline Float64 biquad(Float64 x, Float64 *filt, Float64 *state)
{
    Float64 w = x + filt[3]*state[0] + filt[4]*state[1];
    Float64 y = filt[0]*w + filt[1]*state[0] + filt[2]*state[1];
    state[1] = state[0];
    state[0] = w;
    return y;
}

// -------------------------------------------------------------------------------------
Float64 TFilter::filter(Float64 v)
{
	return biquad(v, m_coffs, m_state);
}

// -------------------------------------------------------------------------------------
TFilter *TFilter::clone()
{
	return new TFilter(m_coffs);
}

// ---------------------------------------------------------------------------
const Float64 PI = acos(-1.0);

void compute_bpf_coeffs(Float64 fsamp, Float64 fc, Float64 q, Float64 gaindB,
			Float64 *pcoffs)
{
  Float64 cf = cos(2.0*PI*fc/fsamp);
  Float64 theta = PI*fc/(fsamp*q);
  Float64 sq = sin(theta);
  Float64 cq = cos(theta);
  Float64 beta = (cq-sq)/(cq+sq);
  Float64 gamma = cf*(1.0+beta);
  Float64 alpha = 0.5*(1.0-beta);
  Float64 ampl = dbampl20(gaindB);
  pcoffs[0] = ampl * alpha;
  pcoffs[1] = 0.0;
  pcoffs[2] = -ampl * alpha;
  pcoffs[3] = gamma;
  pcoffs[4] = -beta;
}

// -------------------------------------------------------------------------------------
void compute_hishelf_coeffs(Float64 fsamp, Float64 fc, Float64 q, Float64 leveldB,
			    Float64 *pcoffs)
{
  Float64 cf = cos(2.0*PI*fc/fsamp);
  Float64 sf = sin(2.0*PI*fc/fsamp);
  Float64 a  = dbampl20(leveldB/2.0);
  Float64 beta = sqrt(a)/q;
  Float64 b0 = a*((a+1.0) + cf*(a-1.0) + sf*beta);
  Float64 b1 = -2.0*a*((a-1.0) + cf*(a+1.0));
  Float64 b2 = a*((a+1.0) + cf*(a-1.0) - sf*beta);
  Float64 a0 = (a+1.0) - cf*(a-1.0) + sf*beta;
  Float64 a1 = 2.0*((a-1.0) - cf*(a+1.0));
  Float64 a2 = (a+1.0) - cf*(a-1.0) - sf*beta;
  pcoffs[0] = b0/a0;
  pcoffs[1] = b1/a0;
  pcoffs[2] = b2/a0;
  pcoffs[3] = -a1/a0;
  pcoffs[4] = -a2/a0;
}

// -------------------------------------------------------------------------------------
void compute_lpf_coeffs(Float64 fsamp, Float64 fc, Float64 q, Float64 gaindB, Float64 *pcoffs)
{
	Float64 w0 = 2.0*PI*fc/fsamp;
	Float64 cf = cos(w0);
	Float64 alpha = sin(w0)/(2.0*q);
	Float64 g  = (1.0 - cf);
	Float64 ampl = dbampl20(gaindB);
	Float64 b0 = ampl * 0.5*g;
	Float64 b1 = ampl * g;
	Float64 b2 = b0;
	Float64 a0 = 1.0 + alpha;
	Float64 a1 = -2.0*cf;
	Float64 a2 = 1.0 - alpha;

	pcoffs[0] = b0/a0;		// b0
	pcoffs[1] = b1/a0;		// b1
	pcoffs[2] = b2/a0;		// b2
	pcoffs[3] = -a1/a0;		// a1
	pcoffs[4] = -a2/a0;		// a2
}

// -------------------------------------------------------------------------------------
void compute_hpf_coeffs(Float64 fsamp, Float64 fc, Float64 q, Float64 gaindB, Float64 *pcoffs)
{
	Float64 w0 = 2.0*PI*fc/fsamp;
	Float64 cf = cos(w0);
	Float64 alpha = sin(w0)/(2.0*q);
	Float64 g  = (1.0 + cf);
	Float64 ampl = dbampl20(gaindB);
	Float64 b0 = ampl * 0.5*g;
	Float64 b1 = -(ampl * g);
	Float64 b2 = b0;
	Float64 a0 = 1.0 + alpha;
	Float64 a1 = -2.0*cf;
	Float64 a2 = 1.0 - alpha;

	pcoffs[0] = b0/a0;		// b0
	pcoffs[1] = b1/a0;		// b1
	pcoffs[2] = b2/a0;		// b2
	pcoffs[3] = -a1/a0;		// a1
	pcoffs[4] = -a2/a0;		// a2
}

// ----------------------------------------------------------------------------

TFilter *make_bpf_filter(Float64 fsamp, Float64 fc, Float64 q, Float64 gaindB)
{
	Float64 coffs[5];
	compute_bpf_coeffs(fsamp, fc, q, gaindB, coffs);
	return new TFilter(coffs);
}

// -------------------------------------------------------------------------------------
TFilter *make_hishelf_filter(Float64 fsamp, Float64 fc, Float64 q, Float64 leveldB)
{
	Float64 coffs[5];
	compute_hishelf_coeffs(fsamp, fc, q, leveldB, coffs);
	return new TFilter(coffs);
}

// -------------------------------------------------------------------------------------
TFilter *make_lpf_filter(Float64 fsamp, Float64 fc, Float64 q, Float64 gaindB)
{
	Float64 coffs[5];
	compute_lpf_coeffs(fsamp, fc, q, gaindB, coffs);
	return new TFilter(coffs);
}

// -------------------------------------------------------------------------------------
TFilter *make_hpf_filter(Float64 fsamp, Float64 fc, Float64 q, Float64 gaindB)
{
	Float64 coffs[5];
	compute_hpf_coeffs(fsamp, fc, q, gaindB, coffs);
	return new TFilter(coffs);
}

// ----------------------------------------------------------

TBWeightedFilter::TBWeightedFilter(Float64 fsamp)
{
    Float64 coffs[5];
    
	Float64 fcx = 1560.0; // best fit from 10 Hz to 10 kHz
	Float64 gamma = tan(fcx * PI/fsamp) / fcx;
	Float64 b0 = 1.0;
	Float64 b1 = -1.0;
	Float64 b2 = 0.0;
	Float64 a0 = 317.0 * gamma + 2.0;
	Float64 a1 = 317.0 * gamma - 2.0;
	Float64 a2 = 0.0;
	Float64 gn = 7.442e9 * gamma * gamma;
    coffs[0] = gn * b0/a0;
    coffs[1] = gn * b1/a0;
    coffs[2] = gn * b2/a0;
    coffs[3] = -a1/a0;
    coffs[4] = -a2/a0;
    m_filter1 = new TFilter(coffs);
    
	b0 = 1.0;
	b1 = 2.0;
	b2 = 1.0;
	a0 = 103.0*gamma + 5.0;
	a2 = 103.0*gamma - 5.0;
	a1 = 2.0*a0*a2;
	a0 *= a0;
	a2 *= a2;
    coffs[0] = b0/a0;
    coffs[1] = b1/a0;
    coffs[2] = b2/a0;
    coffs[3] = -a1/a0;
    coffs[4] = -a2/a0;
    m_filter2 = new TFilter(coffs);
    
	b0 = 1.0;
	b1 = -2.0;
	b2 = 1.0;
	a0 = 12200.0 * gamma + 1.0;
	a2 = 12200.0 * gamma - 1.0;
	a1 = 2.0 * a0 * a2;
	a0 *= a0;
	a2 *= a2;
    coffs[0] = b0/a0;
    coffs[1] = b1/a0;
    coffs[2] = b2/a0;
    coffs[3] = -a1/a0;
    coffs[4] = -a2/a0;
    m_filter3 = new TFilter(coffs);
}

Float64 TBWeightedFilter::filter(Float64 x)
{
    return m_filter3->filter(m_filter2->filter(m_filter1->filter(x)));
}

void TBWeightedFilter::reset()
{
    m_filter1->reset();
    m_filter2->reset();
    m_filter3->reset();
}

#if 1
void TBWeightedFilter::rdz()
{
    m_filter1->rdz();
    m_filter2->rdz();
    m_filter3->rdz();
}
#endif

// ----------------------------------------------------------

TOrd4Filter::TOrd4Filter(Float64 *pcoffs1, Float64 *pcoffs2)
{
    m_filter1 = new TFilter(pcoffs1);
    m_filter2 = new TFilter(pcoffs2);
}

void TOrd4Filter::reset()
{
    m_filter1->reset();
    m_filter2->reset();
}

#if 1
void TOrd4Filter::rdz()
{
    m_filter1->rdz();
    m_filter2->rdz();
}
#endif

Float64 TOrd4Filter::filter(Float64 x)
{
    return m_filter2->filter(m_filter1->filter(x));
}
