// ------------------------------------------------------
// C + Assembly Version
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
#include <float.h>
#ifndef MACOS
#include <time.h>
#endif
#include <stdlib.h>
//#include <dvec.h>
#include "my_types.h"
#include "old-dither.h"

#include "Version.h"

#ifndef MACOS
extern Float32 ran1(SInt32 *idum);
static SInt32  seed = -1;
#endif

#define DITHER_TABLE_SIZE  0x20000	// 128K entries, must be pow(2)

// -------------------------------------------------------------------------------------

#if 0
static const Float32 kDither_double_to_float = ldexpf(1.0f,-23);
static const Float32 kDither_denorm          = ldexpf(1.0f,-34);
#endif

// -------------------------------------------------------------------------------------

TRootDither::TRootDither()
{
#ifndef MACOS
    time((time_t*)&seed);
#endif
    
    fill_dither_table();
}

TRootDither::~TRootDither()
{
    delete dither_table;
}

// ----------------------------------------------------------------------
// TPDF dither

Float32 uniform_variate()
{
    // return a uniform random variate in the range (0.0, 1.0)
#ifdef MACOS
    return ldexpf((Float32)arc4random(),-32);
#else
    return ran1(&seed);
#endif
}

Float32 tpdf_variate()
{
    // return a TPDF random variate in the range (-1.0, 1.0)
    return (uniform_variate() - uniform_variate());
}

void TRootDither::fill_dither_table()
{
    if(!dither_table)
    {
        dither_table_size = DITHER_TABLE_SIZE;
        dither_table = new Float32[dither_table_size];
        
        // Flat spectral shaping
        // the difference of two successive uniform random values
        // produces a flat spectrum, but TPDF values
        // resulting random TPDF has range (-1.0,1.0) ULP.
        
        for(UInt32 ix = 0; ix < dither_table_size; ++ix)
            dither_table[ix] = ldexpf(tpdf_variate(), -24);
        
        dither_index = 0;
    }
}

TRootDither gDither;

// ----------------------------------------------------------------------
TDither::TDither(UInt32 nel)
{
    dither_table = gDither.get_dither_block(nel);
    dither_index = 0;
    dither_table_size = nel;
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// dither_block -- return the address of a block of dither values of length nsamp
// result should be treated as read-only

Float32* TDither::get_dither_block(UInt32 nsamp)
{
    UInt32 offset = dither_index;
    dither_index += nsamp;
    if(dither_index >= dither_table_size)
    {
        dither_index = nsamp;
        offset = 0;
    }
    return (dither_table + offset);
}

// -------------------------------------------------------------------------------------

Float64 TDither::next_dither()
{
    if(dither_table)
    {
        UInt32 ix = dither_index++;
        if(dither_index >= dither_table_size)
            dither_index = 0;
        return (Float64)dither_table[ix];
    }
    else
        return (Float64)ldexpf(tpdf_variate(), -24);
}

// -------------------------------------------------------------------------------------
#if 0
Float64 TDither::dither_dtos()
{
    return next_dither() * kDither_double_to_float;
}

// -------------------------------------------------------------------------------------

Float64 TDither::denorm_dither()
{
    return next_dither() * kDither_denorm;
}
#endif
// ------------------------------------------------
// Non-parallel versions
//

#if 0
#define DITHER_DENORM(pdst,type,psrc,pdith,ix)  \
((pdst)[(ix)] = (type)((psrc)[(ix)] + (pdith)[(ix)] * kDither_denorm))
#else
#define DITHER_DENORM(pdst,type,psrc,pdith,ix)  \
((pdst)[(ix)] = (type)((psrc)[(ix)] + (pdith)[(ix)]))
#endif

void TDither::copy_stod_with_denorm_dither(Float32 *psrc, Float64 *pdst, UInt32 nel)
{
    if(dither_table)
    {
        Float32 *pdith = get_dither_block(nel);
        
        if(1 & nel)
        {
            DITHER_DENORM(pdst, Float64, psrc, pdith, 0);
            ++psrc;
            ++pdst;
            ++pdith;
        }
        if(2 & nel)
        {
            DITHER_DENORM(pdst, Float64, psrc, pdith, 0);
            DITHER_DENORM(pdst, Float64, psrc, pdith, 1);
            psrc  += 2;
            pdst  += 2;
            pdith += 2;
        }
        if(4 & nel)
        {
            DITHER_DENORM(pdst, Float64, psrc, pdith, 0);
            DITHER_DENORM(pdst, Float64, psrc, pdith, 1);
            DITHER_DENORM(pdst, Float64, psrc, pdith, 2);
            DITHER_DENORM(pdst, Float64, psrc, pdith, 3);
            psrc  += 4;
            pdst  += 4;
            pdith += 4;
        }
        nel >>= 3;
        while(nel-- > 0)
        {
            DITHER_DENORM(pdst, Float64, psrc, pdith, 0);
            DITHER_DENORM(pdst, Float64, psrc, pdith, 1);
            DITHER_DENORM(pdst, Float64, psrc, pdith, 2);
            DITHER_DENORM(pdst, Float64, psrc, pdith, 3);
            DITHER_DENORM(pdst, Float64, psrc, pdith, 4);
            DITHER_DENORM(pdst, Float64, psrc, pdith, 5);
            DITHER_DENORM(pdst, Float64, psrc, pdith, 6);
            DITHER_DENORM(pdst, Float64, psrc, pdith, 7);
            psrc  += 8;
            pdst  += 8;
            pdith += 8;
        }
    }
}

void TDither::copy_stos_with_denorm_dither(Float32 *psrc, Float32 *pdst, UInt32 nel)
{
    if(dither_table)
    {
        Float32 *pdith = get_dither_block(nel);
        
        if(1 & nel)
        {
            DITHER_DENORM(pdst, Float32, psrc, pdith, 0);
            ++psrc;
            ++pdst;
            ++pdith;
        }
        if(2 & nel)
        {
            DITHER_DENORM(pdst, Float32, psrc, pdith, 0);
            DITHER_DENORM(pdst, Float32, psrc, pdith, 1);
            psrc  += 2;
            pdst  += 2;
            pdith += 2;
        }
        if(4 & nel)
        {
            DITHER_DENORM(pdst, Float32, psrc, pdith, 0);
            DITHER_DENORM(pdst, Float32, psrc, pdith, 1);
            DITHER_DENORM(pdst, Float32, psrc, pdith, 2);
            DITHER_DENORM(pdst, Float32, psrc, pdith, 3);
            psrc  += 4;
            pdst  += 4;
            pdith += 4;
        }
        nel >>= 3;
        while(nel-- > 0)
        {
            DITHER_DENORM(pdst, Float32, psrc, pdith, 0);
            DITHER_DENORM(pdst, Float32, psrc, pdith, 1);
            DITHER_DENORM(pdst, Float32, psrc, pdith, 2);
            DITHER_DENORM(pdst, Float32, psrc, pdith, 3);
            DITHER_DENORM(pdst, Float32, psrc, pdith, 4);
            DITHER_DENORM(pdst, Float32, psrc, pdith, 5);
            DITHER_DENORM(pdst, Float32, psrc, pdith, 6);
            DITHER_DENORM(pdst, Float32, psrc, pdith, 7);
            psrc  += 8;
            pdst  += 8;
            pdith += 8;
        }
    }
}

// -------------------------------------------------------------------------------------
#if 0
#define DITHER_DTOF(pdst,type,psrc,pdith,ix)  \
((pdst)[(ix)] = (type)((psrc)[(ix)] + (pdith)[(ix)] * kDither_double_to_float))
#else
#define DITHER_DTOF(pdst,type,psrc,pdith,ix)  \
((pdst)[(ix)] = (type)((psrc)[(ix)] + (pdith)[(ix)]))
#endif

void TDither::copy_dtos_with_dither(Float64 *psrc, Float32 *pdst, UInt32 nel)
{
    if(dither_table)
    {
        Float32 *pdith = get_dither_block(nel);
        
        if(1 & nel)
        {
            DITHER_DTOF(pdst,Float32,psrc,pdith,0);
            ++psrc;
            ++pdst;
            ++pdith;
        }
        if(2 & nel)
        {
            DITHER_DTOF(pdst,Float32,psrc,pdith,0);
            DITHER_DTOF(pdst,Float32,psrc,pdith,1);
            psrc  += 2;
            pdst  += 2;
            pdith += 2;
        }
        if(4 & nel)
        {
            DITHER_DTOF(pdst,Float32,psrc,pdith,0);
            DITHER_DTOF(pdst,Float32,psrc,pdith,1);
            DITHER_DTOF(pdst,Float32,psrc,pdith,2);
            DITHER_DTOF(pdst,Float32,psrc,pdith,3);
            psrc  += 4;
            pdst  += 4;
            pdith += 4;
        }
        nel >>= 3;
        while(nel-- > 0)
        {
            DITHER_DTOF(pdst,Float32,psrc,pdith,0);
            DITHER_DTOF(pdst,Float32,psrc,pdith,1);
            DITHER_DTOF(pdst,Float32,psrc,pdith,2);
            DITHER_DTOF(pdst,Float32,psrc,pdith,3);
            DITHER_DTOF(pdst,Float32,psrc,pdith,4);
            DITHER_DTOF(pdst,Float32,psrc,pdith,5);
            DITHER_DTOF(pdst,Float32,psrc,pdith,6);
            DITHER_DTOF(pdst,Float32,psrc,pdith,7);
            psrc  += 8;
            pdst  += 8;
            pdith += 8;
        }
    }
}

// -------------------------------------------------------------------------------------
