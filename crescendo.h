//-------------------------------------------------------------------
// Crescendo.h -- Crescendo for Mac VST and Mac Audio Units
// Stereo plugin which applies Left/Right Crescendo Hearing Corrections
// (C) 2005, Refined Audiometrics Laboratory, LLC. All rights reserved.
// DM/RAL  01/05-11/16
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

#ifndef __CRESCENDO_H__
#define __CRESCENDO_H__

#include "useful_math.h"
#include "ipp_intf.h"
#include "smart_ptr.h"
#include "circbuf.h"
#include "TFilter.h"
#include "hdpheq.h"
#include "vTuningParams.h"

// -------------------------------------------------------------
// The Crescendo 3D Algorithm
// One band of correction for every Bark band
//
#define NFBANDS         25
#define NSUBBANDS		4

// -------------------------------------------------------------
//
struct bark_rec {
    // the following are dynamically updated during processing
    Float64  prev_pwr;
    Float64  mean_pwr;
    Float64  release;
    SInt32   holdctr;
    
    Float64  prev_gain;
    
    // the following set by set_vtuning for each Bark band
    Float64  interp_frac;
    Float64 *pcoff1;
    Float64 *pcoff2;
};

// -------------------------------------------------------------
//
#define SHARED_VAR(type,name) \
type m_##name; \
type get_##name() { return m_##name; } \
void set_##name(type arg) { m_##name = arg; }

// -------------------------------------------------------------
//
class TCrescendo_bark_channel;

class TCrescendo
{
	TPtr<TCrescendo_bark_channel> m_lchan;
	TPtr<TCrescendo_bark_channel> m_rchan;
    
    DZPtr m_DataWindow;
    DZPtr m_HalfWindow;
    TPtr<ipp_fft> m_AudioFFT;
    
    t_EQStruct *m_HdphEQ_basis;
    t_EQStruct *m_PostEQ_basis;
    Float64     m_InvATH[129];
    
    // coalesced EQ's
    Float64*    m_UnifiedEQ;
    Float64*    m_UnifiedEQAmpl;
    
    DZPtr  m_PowerSpectrum;
    DZPtr  m_Filter;
    DZPtr  m_Data;
    
    Float64 m_BarkSpectrum[128];
    Float64 m_BarkGains[128];
    
    float   m_sampleRate;
    float   m_vTuning;
    
    void init_datawin();
    void fill_bark_interpolation_tables();
    void fill_bark_tables();
    void fill_ft_tables();
    void interpolate_eqStruct(t_EQStruct *eqtbl, Float64 *dst, tAmplFn *pfn, bool norm1kHz = true);
    void compute_inverse_ATH_filter();
    void invalidate_unified_filter();
    void ensure_unified_filter();
    
    // --------------------------------------------------------------
    // High level, unified access to FFT cells
    
    inline void get_FT_cell(Float64* ft, UInt32 ix, Float64 &re, Float64 &im)
    {
        m_AudioFFT->get_FT_cell(ft, ix, re, im);
    }
    
    inline void get_FT_DC(Float64* ft, Float64 &re)
    {
        m_AudioFFT->get_FT_DC(ft, re);
    }
    
    inline void get_FT_Nyquist(Float64* ft, Float64 &re)
    {
        m_AudioFFT->get_FT_Nyquist(ft, re);
    }

    
    inline void set_FT_cell(Float64* ft, UInt32 ix, Float64 re, Float64 im)
    {
        m_AudioFFT->set_FT_cell(ft, ix, re, im);
    }
    
    inline void set_FT_DC(Float64* ft, Float64 re)
    {
        m_AudioFFT->set_FT_DC(ft, re);
    }
    
    inline void set_FT_Nyquist(Float64* ft, Float64 re)
    {
        m_AudioFFT->set_FT_Nyquist(ft, re);
    }
    // --------------------------------------------------------------

    
public:
    
    Float64* get_PowerSpectrum()
	{ return m_PowerSpectrum(); }
    
	Float64* get_BarkSpectrum()
	{ return m_BarkSpectrum; }
    
	Float64* get_BarkGains()
	{ return m_BarkGains; }
    
	Float64* get_Data()
	{ return m_Data(); }
    
    SHARED_VAR(UInt32,   blksize);
    SHARED_VAR(UInt32,   hblksize);
    SHARED_VAR(UInt32,   qblksize);
    
	SHARED_VAR(Float64*, HdphEQ);
    SHARED_VAR(Float64*, PostEQ);
    SHARED_VAR(Float64*, PreEQ);
    SHARED_VAR(Float64*, PreEQAmpl);
    
    SHARED_VAR(bool,     Processing);
    // SHARED_VAR(bool,     CorrectionsOnly);
    
    SHARED_VAR(UInt32,   HoldCt);
    SHARED_VAR(Float64,  ReleaseFast);
    SHARED_VAR(Float64,  ReleaseSlow);
    SHARED_VAR(Float64,  GainRelease);
    
    // SHARED_VAR(Float64,  Brightness);
    SHARED_VAR(Float64,  AttendB);
    SHARED_VAR(Float64,  VoldB);
    SHARED_VAR(Float64,  LevelAlpha);
    SHARED_VAR(Float64,  Foldback);
    SHARED_VAR(Float64,  MaxGain);
    SHARED_VAR(Float64,  CaldBSPL);
    SHARED_VAR(Float64,  CaldBFS);
    
    SHARED_VAR(Float64,  selfCalSF);
    
    // ---------------------------------------------
    TCrescendo(Float64 sampleRate);
    virtual ~TCrescendo();
	
    void init(Float64 sampleRate);
    
	void set_headphone(UInt32 ix);
    void set_headphone(t_EQStruct *eqtbl, bool force = false);
    
	void SetSampleRate(Float64 sampleRate);
    
	void set_vtuning(float vtune);
	
    void set_postEQ(UInt32 ix);
    void set_postEQ(t_EQStruct *eqtbl, bool force = false);

    void set_preEQ();

	void render(float *pinL, float *pinR,
                float *poutL, float *poutR,
                UInt32 nel, bool replace,
                tVTuningParams *parms);
    
#ifdef MACOS
	Float64 get_latency();
#else
	UInt32 get_latency();
#endif
    
    Float64 get_power();
	void get_levels(Float64 &lrms, Float64 &rrms);
    
    Float64 convert_dBFS_to_dBSPL(Float64 pdb)
    { return (pdb + m_CaldBSPL - (m_CaldBFS - 3.0)); }
    
	Float64 compute_bark_powers(Float64 *pwr_spectrum, Float64 *bk_pwr);
    void    compute_ft_gains(Float64 *bark_gains, Float64 *ft_buf);
    void    render_samples(Float64 *pin, Float64 *pout, TCrescendo_bark_channel *chan);
    Float64 dbfs_to_dbhl(Float64 pwrfs, Float64 fletch);
    
	void    update_filter(Float64 *pin, TCrescendo_bark_channel *chan);
	void    compute_power_spectrum(Float64 *pin, TCrescendo_bark_channel *chan);
	void    update_bark_powers(TCrescendo_bark_channel *chan);
	void    compute_filter();
	Float64 compute_cumulative_power(Float64 *pwr_spectrum, Float64 *ft_pwr);
    
    void self_calibrate();
};

// -------------------------------------------------------------------------

#define REF_PARENT(type, name) \
type get_##name() { return (type)(m_parent->get_##name()); }

class TCrescendo_bark_channel
{
	TCrescendo *m_parent;
    
	UInt32 m_blksize;
	
	// data for each Bark band
	bark_rec m_bark[NSUBBANDS*NFBANDS+1];
	
	// the buffer used to accumulate input data and scraps
	DZPtr   m_ibuf;
	UInt32  m_ioff;
	UInt32  m_iscrap;
	
	// the output circular buffer used to handle non-pwr-of-2 renderings
	// we always write half frames, but reads can be anything that size
	// or less.
	TPtr<TCircbuf> m_obuf;
    
	Float64   m_level;
    Float64   m_Crest;
    TPtr<TOrd4Filter>  m_crestFilter;
    
public:
	TCrescendo_bark_channel(TCrescendo *parent);
    virtual ~TCrescendo_bark_channel();
	
    int get_ioff()
    { return m_ioff; }
    
    Float64 get_level()
    { return m_level; }
    
    REF_PARENT(UInt32,      blksize);
    REF_PARENT(UInt32,      hblksize);
    REF_PARENT(UInt32,      qblksize);
    
    REF_PARENT(bool,     Processing);
    // REF_PARENT(bool,     CorrectionsOnly);
    
    // REF_PARENT(Float64,  Brightness);
    REF_PARENT(Float64,  AttendB);
    REF_PARENT(Float64,  VoldB);
    REF_PARENT(Float64,  LevelAlpha);
    REF_PARENT(Float64,  Foldback);
    REF_PARENT(Float64,  MaxGain);
    
    REF_PARENT(UInt32,   HoldCt);
    REF_PARENT(Float64,  ReleaseFast);
    REF_PARENT(Float64,  ReleaseSlow);
    REF_PARENT(Float64,  GainRelease);
    
    REF_PARENT(Float64*, BarkSpectrum);
	REF_PARENT(Float64*, BarkGains);
	REF_PARENT(Float64*, Data);
    
    REF_PARENT(Float64, selfCalSF);
    
    void    render_channel(float *pin, float *pout, UInt32 nel, bool replace);
	void    set_vtuning(float vtune);
	void    SetSampleRate(Float64 sampleRate);
	Float64 compute_hcgain(Float64 dbpwr, bark_rec *pbark);
    void    compute_crest_factor(Float64 *pdata, UInt32 nel);
    
	void    compute_bark_gains();
	void    select_data_for_power_estimation(Float64 *pin,
                                             Float64 *pwr_spectrum,
                                             UInt32  hblksize);
	void    select_data_for_filtering(Float64 *pin,
                                      Float64 *data,
                                      UInt32   hblksize);
    
	void update_level(Float64 total_pwr);

    Float64 convert_dBFS_to_dBSPL(Float64 pdb)
    { return m_parent->convert_dBFS_to_dBSPL(pdb); }
    
    void transfer_results_to_output(Float32* pout, UInt32 nel, bool replace);
    
};


// ---------------------------------------------

#endif // __CRESCENDO_H__

// -- end of Crescendo.h -- //
