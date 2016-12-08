// bark_ft.cpp
// DM/RAL  10/07-11/16
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

#include <memory.h>

#include "Crescendo.h"

// ---------------------------------------------------------
// interpolation index and fraction for FFT cell values to Bark channel values
UInt32  ixbk[128+1];
Float64 fxbk[128+1];

// interpolation index and fraction for Bark channel values to FFT cell values
// for the first 129 FFT cells
UInt32  ixft[NSUBBANDS*NFBANDS+3];
Float64 fxft[NSUBBANDS*NFBANDS+3];

// -----------------------------------------------------
// Interpolation routines (Linear)

inline Float64 ft_to_barkd(Float64 *ft_table, UInt32 bark_chan)
{
	UInt32  ix = ixft[bark_chan];
	Float64 fx = fxft[bark_chan];
    return ((1.0 - fx) * ft_table[ix] + fx * ft_table[ix+1]);
}

// -------------------------------------------------------------------------------------
inline Float64 bark_to_ftf(Float64 *bark_table, UInt32 ft_chan)
{
	UInt32  ix = ixbk[ft_chan];
	Float64 fx = fxbk[ft_chan];
    return ((1.0 - fx) * bark_table[ix] + fx * bark_table[ix+1]);
}

// -------------------------------------------------

Float64 TCrescendo::compute_bark_powers(Float64 *pwr_spectrum, Float64 *bk_pwr)
{
    Float64 pwrsum, re, im;
    Float64 ft_pwr[128+1]; // extra one for interpolation routines
	UInt32  ix;
    // Full-scale sinewave should produce FFT amplitudes of 1/2 at +/- freq,
    // for a total power of 1/2 = -3 dB
    // But data windowing will affect the measured peak values.
    // DBM 11/16 - we are now self calibrating - no need for sf
    
    // DC cell has half contribution
    get_FT_DC(pwr_spectrum, re);
    pwrsum = 0.5*re*re*m_UnifiedEQAmpl[0];
    ft_pwr[0] = pwrsum;

    for(ix = 1; ix < 128; ++ix)
    {
        get_FT_cell(pwr_spectrum, ix, re, im);
        pwrsum += (re*re + im*im)*m_UnifiedEQAmpl[ix];
        ft_pwr[ix] = pwrsum;
    }
    
    // just ignore Nyquist contribution
    ft_pwr[128] = pwrsum;

    // At 48 kHz Fsamp, the highest 1/4-Bark bands used are #97 & #98
    //
    // A masking profile of -10 dB/Bark to the low side, and -20 dB/Bark
    // to the high side, integrates over width of 1 bark (-1/2,+1/2) to
    // 0.5 in power. Hence should be approx equiv to summing over 2
    // 1/4-Bark channels -- actually sounds pretty good!
    //
    Float64 ym1 = 0.0;
    Float64 y0  = 0.0;
    for(ix = 0; ix < NSUBBANDS*NFBANDS; ++ix)
    {
        Float64 yp1 = ft_to_barkd(ft_pwr, ix+1);
        bk_pwr[ix] = (yp1 - ym1);
        ym1 = y0;
        y0  = yp1;
    }
    return pwrsum;
}

// -------------------------------------------------------------------------------------
void TCrescendo::compute_ft_gains(Float64 *bark_gains, Float64 *ft_buf)
{
	Float64 ft_gain;
	UInt32  ix;
    
    ft_gain = ampl20(bark_gains[0] + m_UnifiedEQ[0]);
    set_FT_DC(ft_buf, ft_gain);

    for(ix = 1; ix < 128; ++ix)
    {
        ft_gain = bark_to_ftf(bark_gains, ix) + m_UnifiedEQ[ix];
        ft_gain = ampl20(ft_gain);
        set_FT_cell(ft_buf, ix, ft_gain, 0.0);
    }
    
    // zap the frequency zone above audibility
    for(ix = 128; ix < m_hblksize; ++ix)
        set_FT_cell(ft_buf, ix, 0.0, 0.0);

    // just zap the Nyquist contribution
    set_FT_Nyquist(ft_buf, 0.0);
}

// -------------------------------------------------------------------------------------

void TCrescendo::invalidate_unified_filter()
{
    m_UnifiedEQ = 0;
    m_UnifiedEQAmpl = 0;
}

void TCrescendo::ensure_unified_filter()
{
    static Float64 unifiedEQdB[128];
    static Float64 unifiedEQAmpl[128];
    
    if(0 == m_UnifiedEQ)
    {
        Float64 *pATH  = m_InvATH;
        Float64 *pPre  = m_PreEQAmpl;
        Float64 *pHdph = m_HdphEQ;
        
        Float64 *pPreDB  = m_PreEQ;
        Float64 *pPostDB = m_PostEQ;
        
        for(UInt32 ix = 0; ix < 128; ++ix)
        {
            unifiedEQAmpl[ix] = pATH[ix] * pPre[ix] * pHdph[ix];
            unifiedEQdB[ix] = pPreDB[ix] - pPostDB[ix];
        }
        m_UnifiedEQAmpl = unifiedEQAmpl;
        m_UnifiedEQ = unifiedEQdB;
    }
}

// -------------------------------------------------------------------------------------
Float64 cbr(Float64 fkhz)
{
	return (max(0.0, (26.81 / (1.0 + 1.960/fkhz) - 0.53)));
}

// -------------------------------------------------------------------------------------
Float64 inv_cbr(Float64 zbark)
{
	return (1.960 / (26.81 / (zbark + 0.53) - 1.0));
}

// -------------------------------------------------------------------------------------
void TCrescendo::fill_bark_tables()
{
    ixft[0] = 0;
    fxft[0] = 0.0;
	for(int ix = 0; ix < NSUBBANDS*NFBANDS+3; ++ix)
    {
		Float64 zbark = ((Float64)ix)/NSUBBANDS;
		Float64 fkhz = inv_cbr(zbark);
		Float64 cell = fkhz*1.0e3 * m_blksize/m_sampleRate - 0.5;
		if(cell >= 128.0)
		   cell = 127.99;
        else if(cell < 0.0)
            cell = 0;
		UInt32  icell = (SInt32)floor(cell);
		Float64 fcell = cell - icell;
		ixft[ix] = icell;
		fxft[ix] = fcell;
	}
}

// -------------------------------------------------------------------------------------
void TCrescendo::fill_ft_tables()
{
	ixbk[0] = 0;
	fxbk[0] = 0.0;
	for(int ix = 0; ix <= 128; ++ix)
    {
		Float64 fkhz = ix * 1.0e-3 * m_sampleRate / m_blksize;
		Float64 bark = NSUBBANDS * cbr(fkhz);
		UInt32  ibark = (SInt32)floor(bark);
		Float64 fbark = bark - ibark;
		ixbk[ix] = ibark;
		fxbk[ix] = fbark;
    }
}

// -------------------------------------------------------------------------------------
void TCrescendo::fill_bark_interpolation_tables()
{
  fill_bark_tables();
  fill_ft_tables();
}

