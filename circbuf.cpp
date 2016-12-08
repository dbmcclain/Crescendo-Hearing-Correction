// circbuf.cpp -- Circular Buffers
// DM/RAL  06/11

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

#include "circbuf.h"
#include "useful_math.h"
#include "ipp_intf.h"

//-------------------------------------------------------------------
//
TCircbuf::TCircbuf(UInt32 nel)
{
  m_nel = nel;
  m_get = 0;
  m_put = 0;
  m_pdata.allocz(nel);
}

//-------------------------------------------------------------------
//
TCircbuf::~TCircbuf()
{}

//-------------------------------------------------------------------
//
void TCircbuf::put(Float64 *pdata, UInt32 nel)
{
  UInt32 n = min(nel, m_nel - m_put);
  memcpy(m_pdata() + m_put, pdata, n * sizeof(Float64));

  if(nel > n)
    {
      m_put = nel - n;
      memcpy(m_pdata(), pdata + n, m_put*sizeof(Float64));
    }
  else
    incrmod(m_put, n, m_nel);
}

//-------------------------------------------------------------------
//
void TCircbuf::get(Float64 *pdata, UInt32 nel, bool replace)
{
	UInt32 n = min(nel, m_nel - m_get);
	
	if(replace)
		memcpy(pdata, m_pdata() + m_get, n*sizeof(Float64));
	else
		vectorAddD2(m_pdata() + m_get, pdata, n);
	
	if(nel > n)
    {
		m_get = nel - n;
		if(replace)
			memcpy(pdata + n, m_pdata(), m_get*sizeof(Float64));
		else
			vectorAddD2(m_pdata(), pdata + n, m_get);
    }
	else
		incrmod(m_get, n, m_nel);
}

// -- end of circbuf.cpp -- //
