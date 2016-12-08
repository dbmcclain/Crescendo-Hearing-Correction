/*
 *  crescendo_dspPriv.h
 *  crescendo-dsp
 *
 *  Created by David McClain on 7/15/11.
 *  Copyright 2011 Refined Audiometrics Laboratory, LLC. All rights reserved.
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

/* The classes below are not exported */
#pragma GCC visibility push(hidden)

#include "crescendo.h"
#include "crossover.h"

extern void*  RAL_make_headphone_crossover(Float32 sampleRate);
extern void   RAL_discard_headphone_crossover(void *pcross);
extern void   RAL_headphone_crossover_set_sample_rate(void *pcross, Float32 sampleRate);
extern void   RAL_headphone_crossover_process(void* pcross, Float32 *pinL, Float32 *pinR,
                                                  Float32 *poutL, Float32 *poutR, UInt32 nsamp);



// ---------------------------------------------------------------
extern void*  RAL_make_crescendo_processor(Float32 sampleRate);
extern void   RAL_discard_crescendo_processor(void *pcresc);    
extern void   RAL_crescendo_processor_set_sample_rate(void *pcresc, Float32 sampleRate);
extern void   RAL_crescendo_processor_process(void *pcresc, Float32 *pinL, Float32 *pinR,
                                                  Float32 *poutL, Float32 *poutR,
                                                  UInt32 nel, bool replace,
                                                  tVTuningParams *parms);

extern Float64 RAL_crescendo_processor_get_latency(void *pcresc);
extern Float64 RAL_crescendo_processor_get_power(void *pcresc);

// ---------------------------------------------------------------

#pragma GCC visibility pop
