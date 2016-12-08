/*
 *  crossover.cpp
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

#include "useful_math.h"
#include "crossover.h"
#include "old-dither.h"

// -------------------------------------------------------------------------------------
TCrossOver::TCrossOver(Float64 sampleRate)
{
	SetSampleRate(sampleRate);
}

// -------------------------------------------------------------------------------------
TCrossOver::~TCrossOver()
{}

// -------------------------------------------------------------------------------------
void TCrossOver::reset()
{
	m_dlyL->reset();
	m_dlyR->reset();
	m_bpfL->reset();
	m_hishelfL->reset();
	m_bpfR->reset();
	m_hishelfR->reset();
}

// -------------------------------------------------------------------------------------
void TCrossOver::SetSampleRate(Float64 sampleRate)
{
	UInt32 ndly = (UInt32)ceil(0.6*sampleRate/1000.0);
	m_dlyL = new TDelay(ndly);
	m_dlyR = new TDelay(ndly);
	
	m_bpfL     = make_bpf_filter(sampleRate, 5000.0, 1.0, -12.0);
	m_hishelfL = make_hishelf_filter(sampleRate, 750.0, 0.29, -20.0);
	
	m_bpfR     = m_bpfL->clone();
	m_hishelfR = m_hishelfL->clone();

	m_xatten = dbampl20(-6.0);
}

// -------------------------------------------------------------------------------------
void TCrossOver::filter(Float32 *pinL, Float32 *pinR, Float32 *poutL, Float32 *poutR, UInt32 nsamp)
{
    TDither dithL(nsamp);
    TDither dithR(nsamp);
    
	// dereference instance vars into faster locals
	for (int ix = nsamp; --ix >= 0;)
	{
		Float32 fSrcL = *pinL++;
		Float32 fSrcR = *pinR++;
		Float64 fDlyL = m_dlyL->delay(fSrcL);
		Float64 fDlyR = m_dlyR->delay(fSrcR);
        
		Float64 sampLBPF = m_bpfL->filter(fDlyL);
		Float64 sampLHS  = m_hishelfL->filter(fDlyL);
        
		Float64 sampRBPF = m_bpfR->filter(fDlyR);
		Float64 sampRHS  = m_hishelfR->filter(fDlyR);
        
        *poutL++ = dithL.cvt_dtos(fSrcL + m_xatten * (sampRBPF + sampRHS));
        *poutR++ = dithR.cvt_dtos(fSrcR + m_xatten * (sampLBPF + sampLHS));
        
		// for testing
		//poutl[ix] = (float)(0*fSrcL + /* m_xatten * */ (sampLBPF + sampLHS));
		//poutr[ix] = fSrcL; ////
	}
    m_bpfL->rdz();
    m_bpfR->rdz();
    m_hishelfL->rdz();
    m_hishelfR->rdz();
}

