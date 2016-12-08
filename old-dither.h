// dither.h -- interface for dithering routines
// DM/RAL  10/05

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
#ifndef __DITHER_H__
#define __DITHER_H__

#include "my_types.h"

//-------------------------------------------------------------------
class TDither
{
protected:
	Float32 *dither_table;
	UInt32   dither_index;
	UInt32   dither_table_size;
    
	Float64  next_dither();
	Float32 *get_dither_block(UInt32 nel);
    
public:
    TDither() {};
	TDither(UInt32 nel);
    
	Float64 dither_dtos()
    { return next_dither(); }
    
    Float64 denorm_dither()
    { return next_dither(); }
    
	void copy_stod_with_denorm_dither(Float32 *psrc, Float64 *pdst, UInt32 nel);
	void copy_stos_with_denorm_dither(Float32 *psrc, Float32 *pdst, UInt32 nel);
	void copy_dtos_with_dither(Float64 *psrc, Float32 *pdst, UInt32 nel);
    
    Float64 safe_relax(Float64 &accum, Float64 newval, Float64 tc)
    {
        accum += tc * (newval + denorm_dither() - accum);
        return accum;
    }
    
    Float32 cvt_dtos(Float64 x)
    { return (Float32)(x + dither_dtos()); }
};

// -----------------------------------------------------
class TRootDither : public TDither
{
	void fill_dither_table();
    
public:
    TRootDither();
    virtual ~TRootDither();
};

extern TRootDither gDither;

// -----------------------------------------------------

#endif // __DITHER_H__

// -- end of dither.h -- //
