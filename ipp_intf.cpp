// ipp_intf.cpp -- IPP Interface Stuff
// DM/RAL  06/11
// ------------------------------------------------------
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

#include "Version.h"
#include "ipp_intf.h"

//-------------------------------------------------------------------
#ifdef MACOS
// #if ALTIVEC
extern "C"
void *alloc_align16(UInt32 nbytes)
{
    // allocate to 16-byte alignment
    
    // need to allocate enough room so that there will
    // be space for a pointer back to the allocation base address
    // just in front of the actual Float32 array.
    void *p = malloc(nbytes+32);
    if(0 == p) return 0;
    void *q = (void*)(((long)p+19) & ((long)(~15)));
    ((void**)q)[-1] = p;
    return q;
}

void free_align16(void *p)
{
    free(((void**)p)[-1]);
}
// #endif // ALTIVEC
#endif // MACOS

ipp_fft::ipp_fft()
{
	init();
}

ipp_fft::~ipp_fft()
{
	discard();
}

void ipp_fft::init()
{
	m_FFTSpec   = 0;
	m_FFTBuf    = 0;
	m_FFT_Order = 0;
	m_FFT_Flag  = 0;
    m_blkSize   = 0;
#if MACOS
    m_DataBuf   = 0;
    m_hblkSize  = 0;
#endif
}
	
void ipp_fft::discard()
{
	if(m_FFTSpec)
	{
#if WIN32
		ippsFFTFree_R_64f(m_FFTSpec);
		ippsFree(m_FFTBuf);
        
#elif MACOS
        vDSP_destroy_fftsetupD(m_FFTSpec);
        free_align16(m_FFTBuf);
        free_align16(m_DataBuf);
#endif
	}
}

#if MACOS
void ipp_fft::init(UInt32 fft_order, UInt32 flag)
{
	if(m_FFT_Order != fft_order)
    {
        discard();
        m_FFT_Order = fft_order;
        m_blkSize = (1 << fft_order);
        m_hblkSize = (m_blkSize >> 1);
        m_FFTSpec = vDSP_create_fftsetupD(fft_order, FFT_RADIX2);
        m_FFTBuf  = (Float64*)alloc_align16(2*4*m_blkSize*sizeof(Float64));
        m_DataBuf = (Float64*)alloc_align16(m_blkSize*sizeof(Float64));
    }
    m_FFT_Flag = flag;
}

#elif WIN32
void ipp_fft::init(UInt32 fft_order, UInt32 flag)
{
	if(m_FFT_Order != fft_order ||
		m_FFT_Flag != flag)
	{
		discard();
		m_FFT_Order = fft_order;
		m_FFT_Flag  = flag;
		ippsFFTInitAlloc_R_64f(&m_FFTSpec, fft_order, flag, ippAlgHintNone);

		int bufsiz;
		ippsFFTGetBufSize_R_64f(m_FFTSpec, &bufsiz);
		m_FFTBuf = ippsMalloc_8u(bufsiz);

        m_blkSize = (1 << fft_order);
    }
}
#endif

void ipp_fft::fwd(Float64 *buf)
{
#if WIN32
    ippsFFTFwd_RToPack_64f_I(buf, m_FFTSpec, m_FFTBuf); 
    
#elif MACOS
    // points to data for FFT
    DSPDoubleSplitComplex tmp_data = 
    { m_DataBuf, m_DataBuf + m_hblkSize };	
    
	// points to scratch pad for FFT
    DSPDoubleSplitComplex scratch = 
    { m_FFTBuf, m_FFTBuf + 4*m_blkSize };
	
    // points to data for FFT
    DSPDoubleSplitComplex dst = 
    { buf, buf + m_hblkSize };
    
    vDSP_ctozD((DSPDoubleComplex*)buf, 2, &tmp_data, 1, m_hblkSize);
    vDSP_fft_zroptD(m_FFTSpec, &tmp_data, 1, &dst, 1, 
               &scratch, m_FFT_Order, FFT_FORWARD);
    switch(m_FFT_Flag)
    {
        case IPP_FFT_NODIV_BY_ANY:
        case IPP_FFT_DIV_INV_BY_N:
            nspdbMpy1(0.5,buf,m_blkSize);
            break;
            
        case IPP_FFT_DIV_FWD_BY_N:
            nspdbMpy1(0.5/m_blkSize,buf,m_blkSize);
            break;
    }
#endif
}

void ipp_fft::inv(Float64 *buf)
{
#if WIN32
    ippsFFTInv_PackToR_64f_I(buf, m_FFTSpec, m_FFTBuf);
    
#elif MACOS
    DSPDoubleSplitComplex tmp_data = 
    { m_DataBuf, m_DataBuf + m_hblkSize };	
    
	// points to scratch pad for FFT
    DSPDoubleSplitComplex scratch = 
    { m_FFTBuf, m_FFTBuf + 4*m_blkSize };
	
    // points to data for FFT
    DSPDoubleSplitComplex src = 
    { buf, buf + m_hblkSize };
    
    vDSP_fft_zroptD(m_FFTSpec, &src, 1, &tmp_data, 1,
               &scratch, m_FFT_Order, FFT_INVERSE);
    vDSP_ztocD(&tmp_data, 1, (DSPDoubleComplex*)buf, 2, m_hblkSize);
    switch(m_FFT_Flag)
    {
        case IPP_FFT_NODIV_BY_ANY:
        case IPP_FFT_DIV_FWD_BY_N:
            break;
            
        case IPP_FFT_DIV_INV_BY_N:
            nspdbMpy1(1.0/m_blkSize,buf,m_blkSize);
            break;
    }
#endif
}

void ipp_fft::mulSpec(Float64 *src1, Float64 *src2, Float64 *dst)
{
	// assumes src1 is a real-valued spectrum (no imaginary part)
#if WIN32
#if 1
	ippsMulPack_64f(src1, src2, dst, m_blkSize);
#else
    // data stored as R0,R1,I1,R2,I2,...,Rn/2
    dst[0] = src1[0] * src2[0];
    for(int ix = 1; ix < m_blkSize-1; ix += 2)
    {
        Float64 v = src1[ix];
        dst[ix]   = v * src2[ix];
        dst[ix+1] = v * src2[ix+1];
    }
    dst[m_blkSize-1] = 0.0;
#endif
#else // MACOS
	dmul3(src1, src2,            dst,            m_hblkSize);
	dmul3(src1, src2+m_hblkSize, dst+m_hblkSize, m_hblkSize);
#endif
}

#if MACOS
void ippsMulPack_64f(Float64 *src1, Float64 *src2, Float64 *dst, UInt32 nel)
{
	int nel2 = (nel >> 1);
	DSPDoubleSplitComplex csrc1 = {src1+1, src1 + nel2 + 1};
	DSPDoubleSplitComplex csrc2 = {src2+1, src2 + nel2 + 1};
	DSPDoubleSplitComplex cdst  = {dst + 1, dst + nel2 + 1};
	vDSP_zvmulD(&csrc1, 1, &csrc2, 1, &cdst, 1, nel2-1, 1);
	dst[0]    = src1[0] * src2[0];
	dst[nel2] = src1[nel2] * src2[nel2];
}
#endif

// ------------------------------------------------------

// -- end of ipp_intf.cpp --
