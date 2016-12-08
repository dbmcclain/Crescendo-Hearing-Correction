// smart_ptr.h -- Attempt some type safety for C++
// DM/RAL  06/11
// -----------------------------------------

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
#ifndef __SMART_PTR_H__
#define __SMART_PTR_H__

#include <memory.h>
#include "my_types.h"

template <class __T>
class TPtr
{
	__T* m_ptr;
    
public:
	TPtr()
	{ m_ptr = 0; }
    
    TPtr(UInt32 nel)
    { m_ptr = new __T[nel]; }
    
	~TPtr()
	{ discard(); }
    
	void alloc()
	{ discard();
        m_ptr = new __T; }
    
	void alloc(UInt32 nel)
	{ discard();
        m_ptr = new __T[nel]; }
    
	void allocz(UInt32 nel)
	{ discard();
        m_ptr = new __T[nel];
        memset(m_ptr, 0, nel*sizeof(__T)); }
    
	void realloc()
	{ discard();
        m_ptr = new __T; }
    
	void realloc(UInt32 nel)
	{ discard();
        m_ptr = new __T[nel]; }
    
	void reallocz(UInt32 nel)
	{ discard();
        allocz(nel); }
    
	void discard()
	{ 
        if(m_ptr)
        {
            delete m_ptr;
            m_ptr = 0;
        }
    }
    
	__T* operator()()
	{ return m_ptr; }
    
	__T& operator[](UInt32 ix)
	{ return m_ptr[ix]; }
    
	__T* operator->()
	{ return m_ptr; }
    
	__T* operator=(__T* ptr)
    { if(m_ptr != ptr)
        {
            discard();
            m_ptr = ptr;
        }
        return ptr; }
    
	bool operator!()
	{ return (0 == m_ptr); }
};

#endif // __SMART_PTR_H__

// -- end of smart_ptr.h -- //

