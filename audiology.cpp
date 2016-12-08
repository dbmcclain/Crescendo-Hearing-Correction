// audiology.cpp -- Audiology handling for Crescendo
// DM/RAL  06/11
// --------------------------------------------------
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

#include "Crescendo.h"
#define GEN_POLY_TABLES 1
#include "crescendo_polys.h"

//-------------------------------------------------------------------
//
// -------------------------------------------------------------------------------------
inline Float32 clip_to_range(Float32 val, Float32 vmin, Float32 vmax)
{
	return max(vmin, min(val, vmax));
}

// -------------------------------------------------------------------------------------
#if 0
void TCrescendo_bark_channel::set_vtuning(Float32 vtune)
{
    // vtune is dB / Bark
    
	// slope factor is roughly 75 dB over 19 bark, which corresponds to about 6 kHz.
    // Bark 2.5 = 250 Hz. Assume all bands 250 Hz and below need no correction
    
    int start = roundf( 2.5f*NSUBBANDS); // about 250 Hz
    int stop  = roundf(20.0f*NSUBBANDS); // about 6 kHz
    Float32 vtunx = vtune / NSUBBANDS;   // rate per Bark subband
	
    for(int ix = 0; ix <= NSUBBANDS*NFBANDS; ++ix)
	{
		Float32 v = 0.0f;
        if(ix > start)
        {
            if(ix < stop)
                v = vtunx * (ix - start);
            else
                v = vtunx * (stop - start);
        }
		// if the predicted threshold elevation is above our limiting max value
		// then ramp abruptly down to zero to avoid gratuitous noise generation
#if 0
		if(v > 90.0)
			v = 0.0;
		else
#endif
			v = clip_to_range(v, 0.0f, 80.0);

		Float32 df = v / 5.0f;
		UInt32   jx = (UInt32)floor(df);
		df -= jx;

		m_bark[ix].interp_frac = df;
		m_bark[ix].pcoff1 = gfits[jx];
		m_bark[ix].pcoff2 = gfits[jx+1];
    }
}
#else
void TCrescendo_bark_channel::set_vtuning(Float32 vtune)
{
    // vtune is threshold elevation in dBHL at 4kHz
    Float32 slope = 3.575f / NSUBBANDS;
    int start = roundf( 2.5f*NSUBBANDS); // about 250 Hz
    int stop  = roundf(20.0f*NSUBBANDS); // about 6 kHz
	
    for(int ix = 0; ix <= NSUBBANDS*NFBANDS; ++ix)
	{
        int fx = ix;
        if(fx < start)
            fx = start;
        else if(fx > stop)
            fx = stop;
        
        // 17.5 zbark = 4 kHz
        Float32 v = vtune + slope * (fx - 17.5f*NSUBBANDS);
		v = clip_to_range(v, 0.0f, 80.0f);
        
		Float32 df = v / 5.0f;
		UInt32   jx = (UInt32)floor(df);
		df -= jx;
        
		m_bark[ix].interp_frac = df;
		m_bark[ix].pcoff1 = gfits[jx];
		m_bark[ix].pcoff2 = gfits[jx+1];
    }
}
#endif

// -- end of audiology.cpp -- //
