// hdpheq.h
//
// DM/RAL  10/07
// ----------------------------------------

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
#ifndef __HDPHEQ_H__
#define __HDPHEQ_H__


// ------------------------------------- //
// Deemphasis for Allen & Heath GL3300 Console HiShelf EQ
//
// Sample Rate: 48000
// FFT: 256
// Frequency Resolution: 187.5
// Data Window: Hanning
// Y +/-: 0.0
// 256 average(s)
//
//
// ------------------------------------- //
// Deemphasis for Behringer Eurorack UB802
//
// Sample Rate: 48000
// FFT: 256
// Frequency Resolution: 187.5
// Data Window: Hanning
// Y +/-: 0.0
// 256 average(s)
//
//
// extern Float64 gBehrDT880EQ[];

// -----------------------------------------------------------
// --------------- Reference Trace E ---------------
//
// Comment: "Massenburg Average HiShelf 8 kHz"
// Sample Rate: 48000
// FFT: 256
// Frequency Resolution: 187.5
// Data Window: Hanning
// Y +/-: 0.0
// 256 average(s)
//
//

typedef Float64 (tAmplFn)(Float64);

typedef struct t_EQStruct {
	int     fsamp;
	int     nfft;
	Float64 db[129]; // should have NFFT/2+1 elements
} t_EQStruct;

// -------------------------------------------------------------
//

extern t_EQStruct gNullEQ;

#endif // __HDPHEQ_H__
