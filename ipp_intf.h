// ipp_intf.h -- IPP Interface Stuff
// DM/RAL  06/11
// ----------------------------------------------

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
#ifndef __IPP_INTF_H__
#define __IPP_INTF_H__

#include "smart_ptr.h"
#include "memory.h"
#include "Version.h"

// --------------------------------------
#if WIN32xp
#pragma warning(disable : 4305)

extern "C" {
#define nsp_UsesTransform
#define nsp_UsesVector
#include <nsp.h>
};
#endif // WIN32xp

// --------------------------------------
#ifdef MACOS
#include "Accelerate/Accelerate.h"

#define vectorAddD2(src,dst,nel)     vDSP_vaddD(src,1,dst,1,dst,1,nel)
#define copy_dtof(src,dst,nel)       vDSP_vdpsp(src,1,dst,1,nel)
#define copy_ftod(src,dst,nel)       vDSP_vspdp(src,1,dst,1,nel)
#define nspdbMpy3(src1,src2,dst,nel) vDSP_vmulD(src1,1,src2,1,dst,1,nel)
#define nspdbMpy2(src,srcdst,nel)    vDSP_vmulD(src,1,srcdst,1,srcdst,1,nel)
#define nspdbAdd2(src,srcdst,nel)    vDSP_vaddD(src,1,srcdst,1,srcdst,1,nel)
#define nspdbMpy1(k,srcdst,nel)      { Float64 __k = k; vDSP_vsmulD(srcdst,1,&__k,srcdst,1,nel); }

typedef TPtr<Float64> DZPtr;

#if TARGET_CPU_X86
// our compiler does ALL floating point with SSE
#define GETCSR()    ({ int _result; asm volatile ("stmxcsr %0" : "=m" (*&_result) ); /*return*/ _result; })
#define SETCSR( a )    { int _temp = a; asm volatile( "ldmxcsr %0" : : "m" (*&_temp ) ); }

#define DISABLE_DENORMALS int _savemxcsr = GETCSR(); SETCSR(_savemxcsr | 0x8040);
#define RESTORE_DENORMALS SETCSR(_savemxcsr);
#else
#define DISABLE_DENORMALS
#define RESTORE_DENORMALS
#endif

extern void ippsMulPack_64f(Float64 *src1, Float64 *src2, Float64 *dst, UInt32 nel);

#endif // MACOS

// --------------------------------------
#ifdef WIN32

#include <ipp.h>
#include <memory.h>

#define DCplx	  Ipp64fc

#define vectorAddD2(src,dst,nel)  nspdbAdd2(src,dst,nel)
#define nspdbMpy1 ippsMulC_64f_I
#define nspdbMpy2 ippsMul_64f_I
#define nspdbMpy3 ippsMul_64f
#define nspdbAdd2 ippsAdd_64f_I
#define nspsbAdd2 ippsAdd_32f_I
#define nspzbMpy2 ippsMul_64fc_I
#define nspdbCopy ippsCopy_64f
#define nspsbCopy ippsCopy_32f
#define nspsbMpy1 ippsMulC_32f_I
#define nspdMax   ippsMax_64f
#define nspsMax   ippsMax_32f
#define nspdMin   ippsMin_64f
#define nspsMin   ippsMin_32f
#define copy_ftod ippsConvert_32f64f
#define copy_dtof ippsConvert_64f32f

// -----------------------------------------

class DZPtr
{
	Float64* m_ptr;

public:
	DZPtr()
	{ m_ptr = 0; }

	~DZPtr()
	{ ippsFree(m_ptr); }

	void alloc(UInt32 nel)
	{ m_ptr = ippsMalloc_64f(nel); }

	void allocz(UInt32 nel)
	{ m_ptr = ippsMalloc_64f(nel);
	  ippsZero_64f(m_ptr, nel); }

	void realloc(UInt32 nel)
	{ ippsFree(m_ptr);
	   m_ptr = ippsMalloc_64f(nel); }

	void reallocz(UInt32 nel)
	{ ippsFree(m_ptr);
	  allocz(nel); }

	void discard()
	{ ippsFree(m_ptr);
	  m_ptr = 0; }

	Float64* operator()()
	{ return m_ptr; }

	Float64& operator[](UInt32 ix)
	{ return m_ptr[ix]; }
};

// -----------------------------------------
#include <xmmintrin.h>
#include <pmmintrin.h>

class denorm_prot
{
	unsigned int m_dzm;
	unsigned int m_ftz;

public:
	denorm_prot()
	{
		m_ftz = _MM_GET_FLUSH_ZERO_MODE();
		m_dzm = _MM_GET_DENORMALS_ZERO_MODE();
		_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
		_MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
	};

	~denorm_prot()
	{
		_MM_SET_FLUSH_ZERO_MODE(m_ftz);
		_MM_SET_DENORMALS_ZERO_MODE(m_dzm);
	};
};

#endif // WIN32

// -----------------------------------------
// -----------------------------------------

#ifdef MACOS
#define IPP_FFT_DIV_FWD_BY_N   1
#define IPP_FFT_DIV_INV_BY_N   2
#define IPP_FFT_NODIV_BY_ANY   3
#endif

class ipp_fft
{
public:
	ipp_fft();
	virtual ~ipp_fft();
    
	void init(UInt32 fft_order, UInt32 flag = IPP_FFT_DIV_FWD_BY_N);
	void fwd(Float64 *buf);
	void inv(Float64 *buf);
    void mulSpec(Float64 *src1, Float64 *src2, Float64 *dst);

protected:
	virtual void init();
	virtual void discard();
    
	UInt32               m_FFT_Order;
	UInt32               m_FFT_Flag;
    UInt32               m_blkSize;

#if WIN32
	IppsFFTSpec_R_64f	*m_FFTSpec;
	Ipp8u               *m_FFTBuf;

public:
    // --------------------------------------------------------------
    // High level, unified access to FFT cells
    
    inline void get_FT_cell(Float64* ft, int ix, Float64 &re, Float64 &im)
    {
        // warning: index 0 < ix < m_hblksize
        int jx = 2*ix-1;
        re = ft[jx];
        im = ft[jx+1];
    }
    
    inline void get_FT_DC(Float64* ft, Float64 &re)
    {
        re = ft[0];
    }
    
    inline void get_FT_Nyquist(Float64* ft, Float64 &re)
    {
        re = ft[m_blkSize-1]
    }
    
    
    inline void set_FT_cell(Float64* ft, int ix, Float64 re, Float64 im)
    {
        int jx = 2*ix-1;
        ft[jx]   = re;
        ft[jx+1] = im;
    }
    
    inline void set_FT_DC(Float64* ft, Float64 re)
    {
        ft[0] = re;
    }
    
    inline void set_FT_Nyquist(Float64* ft, Float64 re)
    {
        ft[m_blkSize-1] = re;
    }
    // --------------------------------------------------------------

#elif MACOS
    FFTSetupD            m_FFTSpec;
    Float64             *m_DataBuf;
    Float64             *m_FFTBuf;
    UInt32               m_hblkSize;

public:
    // --------------------------------------------------------------
    // High level, unified access to FFT cells
    
    inline void get_FT_cell(Float64* ft, UInt32 ix, Float64 &re, Float64 &im)
    {
        re = ft[ix];
        im = ft[ix+m_hblkSize];
    }
    
    inline void get_FT_DC(Float64* ft, Float64 &re)
    {
        re = ft[0];
    }
    
    inline void get_FT_Nyquist(Float64* ft, Float64 &re)
    {
        re = ft[m_hblkSize];
    }
    
    
    inline void set_FT_cell(Float64* ft, UInt32 ix, Float64 re, Float64 im)
    {
        ft[ix] = re;
        ft[ix+m_hblkSize] = im;
    }
    
    inline void set_FT_DC(Float64* ft, Float64 re)
    {
        ft[0] = re;
    }
    
    inline void set_FT_Nyquist(Float64* ft, Float64 re)
    {
        ft[m_hblkSize] = re;
    }
    // --------------------------------------------------------------
#endif
    
};

inline void dcopy(Float64 *src, Float64 *dst, UInt32 nel)
{
	memcpy(dst, src, nel*sizeof(Float64));
}

inline void dzero(Float64 *p, UInt32 nel)
{
	memset(p, 0, nel*sizeof(Float64));
}

inline void dmul3(Float64 *src1, Float64 *src2, Float64 *dst, UInt32 nel)
{
	nspdbMpy3(src1, src2, dst, nel);
}

inline Float32 dotprodf(Float32 *src1, Float32 *src2, UInt32 nel)
{
	Float32 sum;
#if WIN32
	ippsDotProd_32f(src1, src2, nel, &sum);
#elif MACOS
	vDSP_dotpr(src1, 1, src2, 1, &sum, nel);
#endif
	return sum;
}

#endif // __IPP_INTF_H__

// -- end of ipp_intf.h --
