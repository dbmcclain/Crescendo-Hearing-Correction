/*
 *  tfilter.h
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

#ifndef __TFILTER_H__
#define __TFILTER_H__

#include <math.h>
#include "my_types.h"
#include "smart_ptr.h"

// --------------------------------------------------
class TFilter
{
	Float64 m_coffs[5];
	Float64 m_state[2];
	
public:
	TFilter(Float64 *pcoffs);
    virtual ~TFilter();
	
	Float64 filter(Float64 v);
	void reset();
    void rdz();
    
	TFilter *clone();
};

// --------------------------------------------------
class TOrd4Filter
{
    TPtr<TFilter> m_filter1;
    TPtr<TFilter> m_filter2;
    
public:
    TOrd4Filter(Float64 *pcoffs1, Float64 *pcoffs2);
    
    Float64 filter(Float64);
    void reset();
    void rdz();
};

// --------------------------------------------------
class TBWeightedFilter
{
    TPtr<TFilter> m_filter1;
    TPtr<TFilter> m_filter2;
    TPtr<TFilter> m_filter3;
    
public:
    TBWeightedFilter(Float64 fsamp);
    
    Float64 filter(Float64 x);
    void reset();
    void rdz();
};

// --------------------------------------------------

extern TFilter *make_bpf_filter(Float64 fsamp, Float64 fc, Float64 q, Float64 gaindB);
extern TFilter *make_hishelf_filter(Float64 fsamp, Float64 fc, Float64 q, Float64 leveldB);
extern TFilter *make_lpf_filter(Float64 fsamp, Float64 fc, Float64 q, Float64 gaindB);
extern TFilter *make_hpf_filter(Float64 fsamp, Float64 fc, Float64 q, Float64 gaindB);

inline Float64 dbampl20(Float64 xdb)
{
  return pow(10.0, 0.05*xdb);
}

inline Float32 dbampl20f(Float32 xdb)
{
  return powf(10.0f, 0.05f*xdb);
}

#endif // __TFILTER_H__
