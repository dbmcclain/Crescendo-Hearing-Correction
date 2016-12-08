//-------------------------------------------------------------------
// Crescendo.cpp -- Crescendo for Mac VST and Mac Audio Units
// Stereo plugin which applies Left/Right Crescendo Hearing Corrections
// (C) 2005-2016, Refined Audiometrics Laboratory, LLC. All rights reserved.
// DM/RAL  01/05-12/16
//-------------------------------------------------------------------
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
#include <stdlib.h>
//#include <float.h>

#include "Crescendo.h"
#include "deemph.h"
#include "fletch.h"
#include "crescendo_polys.h"
#include "bark_fft.h"
#include "ipp_intf.h"
#include "old-dither.h"

// -------------------------------------------------------------
//
// -------------------------------------------------------------
//
Float64 e_folding(Float64 tcms, Float64 fsamp)
{
    if(tcms <= 0.0)
        return 1.0;
    else
        return (1.0 - exp(-1.0/(tcms * fsamp)));
}

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------

// ------------------------------------------------------------------
// The Crescendo 3D Algorithm
// FFT-based 1 Correction Band per Critical Band
//

// -------------------------------------------------------------------------------------

void TCrescendo::init_datawin()
{
    // Hann Window over the whole block for power estimation
	Float64 pif = acos(-1.0) / m_blksize;
    
	m_DataWindow.realloc(m_blksize);
    m_DataWindow[0] = 0.0;
	m_DataWindow[m_hblksize] = 1.0;
	for(UInt32 ix = 1; ix < m_hblksize; ++ix)
	{
		Float64 v = sin(pif * ix);
        v *= v;
		m_DataWindow[ix]           = v;
		m_DataWindow[m_blksize-ix] = v;
	}
    
    // Hann Window over half block for filter formation
    m_HalfWindow.realloc(m_hblksize);
    m_HalfWindow[0] = 0.0;
    m_HalfWindow[m_qblksize] = 1.0;
    for(UInt32 ix = 1; ix < m_qblksize; ++ix)
    {
        Float64 v = sin(2.0 * pif * ix);
        v *= v;
        m_HalfWindow[ix]            = v;
        m_HalfWindow[m_hblksize-ix] = v;
    }
}

// ---------------------------------------------------------
Float64 TCrescendo_bark_channel::compute_hcgain(Float64 dbpwr, bark_rec *pbark)
{
	Float64 dbgain = 0.0;
	
	if(dbpwr > 30.0)		// audible?
    {
#if 1
		Float64 foldback = get_Foldback();
		if(dbpwr < foldback)		// range reduction
        {
            dbpwr /= foldback;
            // quadratic ramp to zero from foldback level
            dbpwr = 100.0 - dbpwr * dbpwr * (100.0 - foldback);
        }
#endif
		dbgain = POST_SCALE(INTERPOLATE(dbpwr, pbark->interp_frac,
                                        pbark->pcoff1, pbark->pcoff2));
        dbgain = min(get_MaxGain(), dbgain);
        
#if 1
        // dynamic profile on correction gain
        // to protect against abrupt gain rise around
        // foldback level at noise floor
        
        Float64 gprev = pbark->prev_gain;
        Float64 grls  = get_GainRelease();

        if(dbgain < (gprev - 6.0))
            gprev = dbgain;
        else
            gDither.safe_relax(gprev, dbgain, grls);
        pbark->prev_gain = gprev;
        dbgain = gprev;
#endif
        
    }
	return dbgain;
}

// -------------------------------------------------------------------------------------
Float64 TCrescendo::dbfs_to_dbhl(Float64 pwrfs, Float64 fletch)
{
	// convert from dBFS to dBHL, ATH already removed in FFT power computation
	Float64 pwrhl  = convert_dBFS_to_dBSPL(pwrfs);
	if(fletch < 0.0)
		pwrhl /= 1.0 + fletch / 240.0;
	else if (fletch < 120.0)
		pwrhl /= 1.0 - fletch / 120.0;
	return pwrhl;
}

// -------------------------------------------------------------------------------------
inline Float64 gdbhl_to_gdbspl(Float64 gainhl, Float64 fletch)
{
	// convert from gain in dBHL back to electrical gain in dBSPL
	if(fletch < 0.0)
		gainhl *= 1.0 + fletch / 240.0;
	else
		gainhl *= 1.0 - fletch / 120.0;
	return gainhl;
}

// -------------------------------------------------------------------------------------

void TCrescendo::update_filter(Float64 *pin, TCrescendo_bark_channel *chan)
{
    compute_power_spectrum(pin, chan);
    update_bark_powers(chan);
    chan->compute_bark_gains();
    compute_filter();
}

void TCrescendo_bark_channel::compute_crest_factor(Float64 *pdata, UInt32 nel)
{
    Float64 pkval = 0.0;
    Float64 sumsq = 0.0;
    for(int ix = 0; ix < nel; ++ix)
    {
        Float64 v = pdata[ix];
        v = m_crestFilter->filter(v);
        v *= v;
        sumsq += v;
        if(v > pkval)
            pkval = v;
    }
    // TC is 200 ms at 48 kHz
    m_Crest = db10(pkval) - db10(sumsq/nel);
    m_crestFilter->rdz();
}


void TCrescendo_bark_channel::select_data_for_power_estimation(Float64 *pin,
                                                               Float64 *pwr_spectrum,
                                                               UInt32  hblksize)
{
    // windowed transform of latest half-block for filter creation
    switch(m_ioff)
    {
        case 0:
            dcopy(pin+2*hblksize, pwr_spectrum,          hblksize);
            dcopy(pin+0*hblksize, pwr_spectrum+hblksize, hblksize);
            break;
            
            
        case 1:
            dcopy(pin+0*hblksize, pwr_spectrum,          2*hblksize);
            break;
            
        case 2:
            dcopy(pin+1*hblksize, pwr_spectrum,          2*hblksize);
            break;
    }
}

void TCrescendo::compute_power_spectrum(Float64 *pin, TCrescendo_bark_channel *chan)
{
    Float64 *pwr_spectrum = get_PowerSpectrum();
    Float64 *pwin = m_DataWindow();
    
    chan->select_data_for_power_estimation(pin, pwr_spectrum, m_hblksize);
    chan->compute_crest_factor(pwr_spectrum + m_qblksize, m_hblksize);
    
    dmul3(pwin, pwr_spectrum, pwr_spectrum, (int)m_blksize);
    m_AudioFFT->fwd(pwr_spectrum);
}

void TCrescendo::self_calibrate()
{
    // Self calibration for power estimation
    // The character of the FFT routine, its scaling,
    // and the effects of data windowing will modify the
    // result from expected ideal values.
    //
    // This routine measures that variation and generates a
    // correction term for use during the running measurements
    // of signal power.
    
    Float64 *pwr_spectrum = get_PowerSpectrum();
    Float64 *pwin = m_DataWindow();
    Float64 *bark_spectrum = get_BarkSpectrum();
    
    // Construct a unit amplitude sinewave exactly subharmonic of sample rate.
    // At 48 kHz, this is about 1.5 kHz.

    Float64 pif = 2.0*acos(-1.0)*32.0/m_blksize;
    
    for(UInt32 ix = 0; ix < m_blksize; ++ix)
        pwr_spectrum[ix] = sin(pif * ix);
    
    // build a dummy EQ amplitude table for unity gain processing
    Float64 dummyEQ[128];
    for(int ix = 0; ix < 128; ++ix)
        dummyEQ[ix] = 1.0;
    
    // Sum over the signal power the same way we do during processing.

    dmul3(pwin, pwr_spectrum, pwr_spectrum, (int)m_blksize);
    m_AudioFFT->fwd(pwr_spectrum);
    
    Float64 *savEQ = m_UnifiedEQAmpl;
    m_UnifiedEQAmpl = dummyEQ;
    Float64 pwrsum = compute_bark_powers(pwr_spectrum, bark_spectrum);
    m_UnifiedEQAmpl = savEQ;
    
    // save the factor which represents twice this total power for the 0 dB level.
    m_selfCalSF = db10(2.0*pwrsum);
}


void TCrescendo::update_bark_powers(TCrescendo_bark_channel *chan)
{
    Float64 *bark_spectrum = get_BarkSpectrum();
    Float64 *pwr_spectrum  = get_PowerSpectrum();
    
    // split FFT power into Bark bands
    Float64 total_pwr = db10(compute_bark_powers(pwr_spectrum, bark_spectrum)) - get_selfCalSF();
    
    // compute ongoing level for calibration (300 ms)
    chan->update_level(total_pwr);
}

void TCrescendo_bark_channel::update_level(Float64 total_pwr)
{
    gDither.safe_relax(m_level, total_pwr, get_LevelAlpha());
}

void TCrescendo_bark_channel::compute_bark_gains()
{
	Float64 *bpwr  = get_BarkSpectrum();
	Float64 *bgain = get_BarkGains();
    
	// compute upward and downward masking contributions
	//compute_masks(cpwr);
    
	// now compute Bark channel output gains
	bark_rec *pbark     = m_bark;
	Float64 attendb     = get_AttendB(); // processing headroom, to be made up externally
	Float64 voldb       = get_VoldB();   // desired volume boost
    Float64 gain0       = attendb + voldb;
	Float64 releaseSlow = get_ReleaseSlow();
	Float64 releaseFast = get_ReleaseFast();
	UInt32  holdct      = get_HoldCt();
    
	// Float64 brightness = get_brightness();
	bool   proc        = get_Processing();
    // bool   corr_only   = get_CorrectionsOnly();
    
#if 0
	// 250 Hz = Bark 2.5 = 1/4-Bark band #10
	// so < 259 Hz we have no VTuning corrections
	for(int ix = 0; ix < 10; ++ix, ++pbark)
	{
		// bgain[ix] = attendb + pbark->pre_eq - pbark->hdph + voldb;
		bgain[ix] = corr_only ? -200.0 : gain0;
	}
#endif
    // > 250 Hz we have VTuning corrections
	for(int ix = 0; ix < NSUBBANDS*NFBANDS; ++ix, ++pbark)
	{
		// -----------------------------------------------------------------------
		// compute attack and release on measured power
		// uses a fast attack and slower release, plus a hold
        
		Float64 xpwr = db10(bpwr[ix]) - get_selfCalSF();
		Float64 mn   = pbark->mean_pwr;
		gDither.safe_relax(mn, xpwr, releaseSlow);
		pbark->mean_pwr = mn;
        
		Float64 crest = m_Crest;
		Float64 prev = pbark->prev_pwr;
        
		if(xpwr > prev)
		{
            // instantaneous attack >= 6 dB
            if(xpwr + crest > prev + 6.0)
		    {
                pbark->release = releaseFast; // 50 ms
                pbark->holdctr = holdct;
                prev = xpwr + crest;
		    }
            else
                gDither.safe_relax(prev, xpwr, releaseFast);
		}
		else if(pbark->holdctr > 0)
		{
			// hold until hold counter empties
			--pbark->holdctr;
		}
		else
		{
			// if releasing from strong impulse (> 6 dB above mean levels)
			// then use a fast release until we fall below 1 dB above mean levels
			//
			// Otherwise if we are within 3 dB of mean levels use a slow release
			// We have hysteresis to avoid toggling between fast and slow release.
			//
			// If we are in-between these two cases, then we continue using the release
			// we last used.
			//
            
			if (prev < mn + 3.0) // 3 dB
				pbark->release = releaseSlow; // 200 ms
            
			gDither.safe_relax(prev, xpwr, pbark->release);
		}
		pbark->prev_pwr = prev;
		xpwr = prev;
        
		// -----------------------------------------------------------------------
		// now compute HC gains...
        
		// the dB power in the cochlea, including masking contributions
		// and fast attack, slow release with hold
        Float64 dgain = 0.0;
        
		if(proc)
		{
			// amount of HC gain needed
            
			// power in the music itself for limiter use
			Float64 pwrdb = xpwr + voldb;
            
			// fletch is just the audible threshold measured relative to 0 dBSPL at 1 kHz
			Float64 fletch = gFletch[ix];
            
			// -----------------------------------------------------------------------
			// conversion from dBFS to dBSPL
            
			Float64 pwrhl = m_parent->dbfs_to_dbhl(pwrdb, fletch);
			Float64 dbhc = compute_hcgain(pwrhl, pbark);
            
			// correction gain to be applied in SPL space
			dgain = gdbhl_to_gdbspl(dbhc, fletch);
		}
#if 0
        if(corr_only)
            bgain[ix] = gain0 + db20(ampl20(dgain) - 1.0);
        else
            bgain[ix] = gain0 + dgain;
#else
        bgain[ix] = gain0 + dgain;
#endif
	}
}

void TCrescendo::compute_filter()
{
    Float64 *bgain = get_BarkGains();
    Float64 *pwr_spectrum = get_PowerSpectrum();
    Float64 *filter = m_Filter();
    Float64 *pwin   = m_HalfWindow();
    
    compute_ft_gains(bgain, pwr_spectrum);
    m_AudioFFT->inv(pwr_spectrum);
    dmul3(pwin+m_qblksize, pwr_spectrum, filter, m_qblksize);
    dzero(filter+m_qblksize, m_hblksize);
    dmul3(pwin, pwr_spectrum+3*m_qblksize, filter+3*m_qblksize, m_qblksize);
    m_AudioFFT->fwd(filter);
}

//-------------------------------------------------------------------

void TCrescendo_bark_channel::select_data_for_filtering(Float64 *pin,
                                                        Float64 *data,
                                                        UInt32   hblksize)
{
    // pin points to a buffer of 3*hblksize input samples
    // m_ioff is the index of half blocks into pin, at which the latest input had been saved
    // delay is 1/2 block = 2.67 ms from input @ 48 kHz
    switch(m_ioff)
    {
        case 0:
            // we just saved at index 0, so process segments 1 & 2
            dcopy(pin+1*hblksize, data,          2*hblksize);
            break;
            
            
        case 1:
            // we just saved at index 1, so process segments 2 & 0
            dcopy(pin+2*hblksize, data,          hblksize);
            dcopy(pin+0*hblksize, data+hblksize, hblksize);
            break;
            
        case 2:
            // we just saved at index 2, so process segments 0 & 1
            dcopy(pin+0*hblksize, data,          2*hblksize);
            break;
    }
}

void TCrescendo::render_samples(Float64 *pin, Float64 *data, TCrescendo_bark_channel *chan)
{
    update_filter(pin, chan);
    
    // non-windowed transform for data
    // overlap-save convolution does not use data windowing
    // 1/2 block delay from filter center = 2.67 ms at 48 kHz
    chan->select_data_for_filtering(pin, data, m_hblksize);
    m_AudioFFT->fwd(data);
    m_AudioFFT->mulSpec(m_Filter(), data, data);
    m_AudioFFT->inv(data);
}

// -------------------------------------------------------------------------------------

void TCrescendo_bark_channel::transfer_results_to_output(Float32* pout, UInt32 nel, bool replace)
{
    Float64 *data = get_Data();
    TDither dith(nel);
    
    m_obuf->get(data, nel, true);
    if(replace)
        dith.copy_dtos_with_dither(data, pout, nel);
    else
    {
        for(UInt32 ix = 0; ix < nel; ++ix)
            pout[ix] = dith.cvt_dtos(pout[ix] + data[ix]);
    }

}

void TCrescendo_bark_channel::render_channel(Float32 *pin,
                                             Float32 *pout,
                                             UInt32   nsamp,
                                             bool     replace)
{
#if 0
	FILE *fp = fopen("/Users/david/Crescendo/MacAudio/cresc_log.txt", "a");
	fprintf(fp, "nsamp = %d\n", nsamp);
	fprintf(fp, "pin = %p\n", pin);
	fprintf(fp, "pout = %p\n", pout);
	fclose(fp);
#endif
	Float64 *data = get_Data();
	int      hblksize = get_hblksize();
    int      qblksize = get_qblksize();
	UInt32   nel  = hblksize - m_iscrap;
    
	while(nsamp >= nel)
    {
        copy_ftod(pin, m_ibuf()+m_ioff*hblksize+m_iscrap, nel);
        
        // we are the half block filled starting at the half-block index m_ioff
		m_parent->render_samples(m_ibuf(), data, this); // results in data
		m_obuf->put(data+qblksize, hblksize);
        transfer_results_to_output(pout, nel, replace);
		
		pin   += nel;
		pout  += nel;
		nsamp -= nel;
		m_iscrap = 0;
		nel = hblksize;
		incrmod(m_ioff, 1, 3);
    }
	if(nsamp > 0)
    {
		copy_ftod(pin, m_ibuf()+m_ioff*hblksize+m_iscrap, nsamp);
		m_iscrap += nsamp;
        transfer_results_to_output(pout, nsamp, replace);
    }
}

// -------------------------------------------------------------------------------------

TCrescendo_bark_channel::TCrescendo_bark_channel(TCrescendo *parent)
{
    int ix;
#if 0
    // ----------------------------------
    fprintf(stderr,"Initialize for sample rate: %F\n", sampleRate);
    fflush(stderr);
    // ----------------------------------
#endif
    
    m_parent = parent;
    m_blksize = 0;
    
    for(ix = NSUBBANDS*NFBANDS+1; --ix >= 0;)
    {
        m_bark[ix].holdctr  = 0;
        m_bark[ix].prev_pwr = -140.0;
        m_bark[ix].mean_pwr = -140.0;
        m_bark[ix].release  = get_ReleaseFast();
        m_bark[ix].prev_gain = 0.0;
    }
    
    set_vtuning(0.0);
    
    // m_limit is max permissible output amplitude
    // Since worst case is 2 consecutive large, equal amplitude, spikes
    // with center amplitude at 4/pi * amplitude of spikes,
    // we should limit max amplitude to less than pi/4 = 0.785 = -2.1 dBFS
    //
    // UPDATE: Research into DAC behavior shows that worst case overshoot/clipping
    // will occur when all the FIR taps abs values add up, in the output
    // interpolation filters. Suggested range of max overshoot is 2.9-5.7dB.
    //
    // Since we can't know how good or bad our filters are,
    // let's try to keep things below -6 dBFS = 0.5.
    
    m_level = 0.0;
    
    // 4th order 500 Hz Butterworth HPF, 48 kHz Fsamp
    static Float64 cfilt_coffs1[5] = {
        0.974538047066236,
        -1.949076094132472,
        0.974538047066236,
        1.9469872972295479,
        -0.9511648910353966};
    static Float64 cfilt_coffs2[5] = {
        0.9420089366880302,
        -1.8840178733760605,
        0.9420089366880302,
        1.8819987984354725,
        -0.8860369483166486};
    m_crestFilter = new TOrd4Filter(cfilt_coffs1, cfilt_coffs2);
}

TCrescendo_bark_channel::~TCrescendo_bark_channel()
{}

//-------------------------------------------------------------------
//
void TCrescendo_bark_channel::SetSampleRate(Float64 sampleRate)
{
#if 0
	// ----------------------------------
	fprintf(stderr,"Setting sample rate: %F\n", sampleRate);
	fflush(stderr);
	// ----------------------------------
#endif
	if(m_parent->m_blksize != m_blksize)
	{
		m_blksize = m_parent->m_blksize;
		int hblksize = m_parent->m_hblksize;
        
		m_ioff   = 0;
		m_iscrap = 0;
        m_ibuf.reallocz(3*hblksize);
		
        m_obuf = new TCircbuf(2*hblksize);
		m_obuf->put(m_ibuf(), hblksize);
    }
}

//-------------------------------------------------------------------
//
TCrescendo::TCrescendo(Float64 sampleRate)
{
    // m_Brightness = 0.0;
	m_AttendB    = 0.0;
	m_VoldB      = 0.0;
	
    m_Foldback = 20.0; // 30.0;
    m_MaxGain  = 50.0; // 40.0; // 24.0;
    m_CaldBFS  = -17.0;
    m_CaldBSPL = 77.0;
    
	m_Processing = true;
	// m_CorrectionsOnly = false;
	
    m_blksize = 0;
    m_sampleRate = 0.0;
    m_vTuning = 0.0;
    
	m_lchan = new TCrescendo_bark_channel(this);
	m_rchan = new TCrescendo_bark_channel(this);
    
    set_vtuning(20.0);
    m_HdphEQ_basis = &gNullEQ;
    m_PostEQ_basis = &gNullEQ;
    invalidate_unified_filter();
    
	SetSampleRate(sampleRate);
}

TCrescendo::~TCrescendo()
{}

//-------------------------------------------------------------------
//
void TCrescendo::SetSampleRate(Float64 sampleRate)
{
    if(m_sampleRate != sampleRate)
    {
        int blksize = (sampleRate > 50.0e3) ? 512 : 256;
        m_sampleRate = sampleRate;
        
        if(blksize != m_blksize)
        {
            m_blksize  = blksize;
            m_hblksize = blksize >> 1;
            m_qblksize = blksize >> 2;
            
            int fft_order = (sampleRate > 50.0e3) ?   9 :   8;
            m_AudioFFT.realloc();
            m_AudioFFT->init(fft_order, IPP_FFT_DIV_FWD_BY_N);
            
            init_datawin();
            
            m_PowerSpectrum.realloc(blksize);
            m_Filter.realloc(blksize);
            m_Data.realloc(blksize);
        }
        
        Float64 blk_sr = sampleRate / m_hblksize;
        m_HoldCt       = (UInt32)ceil(10.0e-3 * blk_sr);
        m_ReleaseSlow  = e_folding(200.0e-3, blk_sr);
        m_ReleaseFast  = e_folding( 50.0e-3, blk_sr);
        m_LevelAlpha   = e_folding(300.0e-3, blk_sr);
        m_GainRelease  = e_folding( 10.0e-3, blk_sr);
        
        fill_bark_interpolation_tables();
        
        compute_inverse_ATH_filter();
        set_preEQ();
        set_postEQ(m_PostEQ_basis, true);
        set_headphone(m_HdphEQ_basis, true);
        ensure_unified_filter();
        
        m_lchan->SetSampleRate(sampleRate);
        m_rchan->SetSampleRate(sampleRate);
        
        self_calibrate();
    }
}

void TCrescendo::set_vtuning(float vtune)
{
    if(m_vTuning != vtune)
    {
        m_vTuning = vtune;
        m_lchan->set_vtuning(vtune);
        m_rchan->set_vtuning(vtune);
    }
}

void TCrescendo::render(Float32 *pinL, Float32 *pinR,
                        Float32 *poutL, Float32 *poutR,
                        UInt32 nel, bool replace,
                        tVTuningParams *parms)
{
    DAZFZ env;
    
    if(parms)
    {
        set_CaldBFS(parms->CaldBFS);
        set_CaldBSPL(parms->CaldBSPL);
        set_vtuning(parms->vTune);
        set_Processing(parms->proc_onoff);
        // set_CorrectionsOnly(parms->corr_only_onoff);
        set_VoldB(parms->voldB);
        set_AttendB(parms->attendB);
        // set_Brightness(parms->brighnessdB);
        set_postEQ(parms->postEQ);
        set_headphone(parms->headphone);
        
        ensure_unified_filter();
    }
	if(pinL && poutL)
		m_lchan->render_channel(pinL, poutL, nel, replace);
	if(pinR && poutR && (pinL != pinR) && (poutL != poutR))
		m_rchan->render_channel(pinR, poutR, nel, replace);
}

void TCrescendo::get_levels(Float64 &lrms, Float64 &rrms)
{
	lrms = m_lchan->get_level();
	rrms = m_rchan->get_level();
}

Float64 TCrescendo::get_power()
{ return m_lchan->get_level(); }

#ifdef MACOS
Float64 TCrescendo::get_latency()
{ return (((Float64)(5*m_qblksize))/m_sampleRate); }
#else
UInt32 TCrescendo::get_latency()
{ return (5*m_qblksize); }
#endif

// -- end of Crescendo.cpp -- //
