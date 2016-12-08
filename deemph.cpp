// deemph.cpp
// DM/RAL  10/07
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
#include "hdpheq.h"

// -------------------------------------------------------------
//
// --------------- Reference Trace A ---------------
//
// Comment: "A&H HiShelf EQ Full Up"
// Sample Rate: 48000
// FFT: 256
// Frequency Resolution: 187.5
// Data Window: Hanning
// Y +/-: 0.0
// 16 average(s)
//
//
// Frequency (Hz)	Magnitude (dB)	Phase (degrees)

#define AHEQ(freq,ampl,phase) (ampl)

t_EQStruct gAH3300EQ = {		// first 129 FFT Cells
    48000, 256,
    {
        AHEQ(0.0,	0.63,	5.85),
        AHEQ(187.5,	0.62,	8.50),
        AHEQ(375.0,	0.72,	11.07),
        AHEQ(562.5,	1.01,	14.08),
        AHEQ(750.0,	1.14,	15.65),
        AHEQ(937.5,	1.43,	19.15),
        AHEQ(1125.0,	1.85,	22.08),
        AHEQ(1312.5,	2.24,	24.22),
        AHEQ(1500.0,	2.52,	25.34),
        AHEQ(1687.5,	2.89,	26.96),
        AHEQ(1875.0,	3.23,	28.49),
        AHEQ(2062.5,	3.50,	29.62),
        AHEQ(2250.0,	3.78,	30.60),
        AHEQ(2437.5,	4.03,	31.51),
        AHEQ(2625.0,	4.27,	32.45),
        AHEQ(2812.5,	4.53,	33.22),
        AHEQ(3000.0,	4.75,	33.95),
        AHEQ(3187.5,	5.01,	34.83),
        AHEQ(3375.0,	5.28,	35.51),
        AHEQ(3562.5,	5.52,	35.95),
        AHEQ(3750.0,	5.75,	36.55),
        AHEQ(3937.5,	5.98,	37.13),
        AHEQ(4125.0,	6.18,	37.54),
        AHEQ(4312.5,	6.35,	37.93),
        AHEQ(4500.0,	6.56,	38.37),
        AHEQ(4687.5,	6.77,	38.91),
        AHEQ(4875.0,	6.99,	39.30),
        AHEQ(5062.5,	7.19,	39.60),
        AHEQ(5250.0,	7.39,	39.90),
        AHEQ(5437.5,	7.57,	40.09),
        AHEQ(5625.0,	7.74,	40.02),
        AHEQ(5812.5,	7.91,	40.10),
        AHEQ(6000.0,	8.07,	40.23),
        AHEQ(6187.5,	8.25,	40.28),
        AHEQ(6375.0,	8.42,	40.37),
        AHEQ(6562.5,	8.59,	40.43),
        AHEQ(6750.0,	8.73,	40.40),
        AHEQ(6937.5,	8.89,	40.38),
        AHEQ(7125.0,	9.05,	40.37),
        AHEQ(7312.5,	9.21,	40.30),
        AHEQ(7500.0,	9.36,	40.20),
        AHEQ(7687.5,	9.52,	40.11),
        AHEQ(7875.0,	9.66,	40.03),
        AHEQ(8062.5,	9.80,	39.82),
        AHEQ(8250.0,	9.93,	39.65),
        AHEQ(8437.5,	10.05,	39.47),
        AHEQ(8625.0,	10.18,	39.31),
        AHEQ(8812.5,	10.31,	39.11),
        AHEQ(9000.0,	10.44,	39.01),
        AHEQ(9187.5,	10.56,	38.82),
        AHEQ(9375.0,	10.69,	38.63),
        AHEQ(9562.5,	10.80,	38.45),
        AHEQ(9750.0,	10.91,	38.22),
        AHEQ(9937.5,	11.01,	37.96),
        AHEQ(10125.0,	11.11,	37.71),
        AHEQ(10312.5,	11.22,	37.42),
        AHEQ(10500.0,	11.32,	37.17),
        AHEQ(10687.5,	11.41,	36.93),
        AHEQ(10875.0,	11.51,	36.69),
        AHEQ(11062.5,	11.61,	36.49),
        AHEQ(11250.0,	11.69,	36.27),
        AHEQ(11437.5,	11.77,	35.93),
        AHEQ(11625.0,	11.86,	35.62),
        AHEQ(11812.5,	11.95,	35.35),
        AHEQ(12000.0,	12.02,	35.07),
        AHEQ(12187.5,	12.12,	34.80),
        AHEQ(12375.0,	12.22,	34.55),
        AHEQ(12562.5,	12.30,	34.36),
        AHEQ(12750.0,	12.38,	34.07),
        AHEQ(12937.5,	12.46,	33.75),
        AHEQ(13125.0,	12.53,	33.46),
        AHEQ(13312.5,	12.59,	33.20),
        AHEQ(13500.0,	12.67,	32.80),
        AHEQ(13687.5,	12.74,	32.47),
        AHEQ(13875.0,	12.81,	32.20),
        AHEQ(14062.5,	12.87,	31.89),
        AHEQ(14250.0,	12.94,	31.60),
        AHEQ(14437.5,	13.01,	31.33),
        AHEQ(14625.0,	13.06,	31.05),
        AHEQ(14812.5,	13.13,	30.82),
        AHEQ(15000.0,	13.20,	30.52),
        AHEQ(15187.5,	13.27,	30.27),
        AHEQ(15375.0,	13.32,	29.97),
        AHEQ(15562.5,	13.38,	29.67),
        AHEQ(15750.0,	13.43,	29.27),
        AHEQ(15937.5,	13.48,	28.95),
        AHEQ(16125.0,	13.52,	28.59),
        AHEQ(16312.5,	13.57,	28.30),
        AHEQ(16500.0,	13.62,	28.00),
        AHEQ(16687.5,	13.67,	27.73),
        AHEQ(16875.0,	13.72,	27.46),
        AHEQ(17062.5,	13.77,	27.19),
        AHEQ(17250.0,	13.82,	26.87),
        AHEQ(17437.5,	13.87,	26.57),
        AHEQ(17625.0,	13.90,	26.27),
        AHEQ(17812.5,	13.95,	25.96),
        AHEQ(18000.0,	13.99,	25.63),
        AHEQ(18187.5,	14.03,	25.30),
        AHEQ(18375.0,	14.06,	24.94),
        AHEQ(18562.5,	14.10,	24.68),
        AHEQ(18750.0,	14.14,	24.43),
        AHEQ(18937.5,	14.18,	24.19),
        AHEQ(19125.0,	14.22,	23.93),
        AHEQ(19312.5,	14.27,	23.70),
        AHEQ(19500.0,	14.30,	23.40),
        AHEQ(19687.5,	14.34,	23.05),
        AHEQ(19875.0,	14.38,	22.73),
        AHEQ(20062.5,	14.42,	22.46),
        AHEQ(20250.0,	14.45,	22.20),
        AHEQ(20437.5,	14.49,	21.86),
        AHEQ(20625.0,	14.52,	21.62),
        AHEQ(20812.5,	14.55,	21.32),
        AHEQ(21000.0,	14.58,	21.00),
        AHEQ(21187.5,	14.61,	20.71),
        AHEQ(21375.0,	14.64,	20.52),
        AHEQ(21562.5,	14.68,	20.25),
        AHEQ(21750.0,	14.71,	19.99),
        AHEQ(21937.5,	14.73,	19.77),
        AHEQ(22125.0,	14.79,	19.54),
        AHEQ(22312.5,	14.86,	19.36),
        AHEQ(22500.0,	14.94,	18.92),
        AHEQ(22687.5,	14.99,	17.77),
        AHEQ(22875.0,	15.01,	14.62),
        AHEQ(23062.5,	14.87,	13.25),
        AHEQ(23250.0,	14.50,	11.30),
        AHEQ(23437.5,	14.45,	8.70),
        AHEQ(23625.0,	14.58,	5.98),
        AHEQ(23812.5,	14.51,	6.61),
        AHEQ(24000.0,	14.61,	4.71) }};

// -------------------------------------------------------------
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
// Frequency (Hz)	Magnitude (dB)	Phase (degrees)
#define behreq(freq, ampdb, phdeg)  (ampdb)
t_EQStruct gBehringerEQ = {
    48000, 256,
    {
        behreq(0.0,	-1.15,	7.83),
        behreq(187.5,	-1.25,	12.24),
        behreq(375.0,	-1.09,	17.81),
        behreq(562.5,	-0.84,	22.17),
        behreq(750.0,	-0.55,	26.07),
        behreq(937.5,	0.01,	33.73),
        behreq(1125.0,	0.75,	40.06),
        behreq(1312.5,	1.48,	44.99),
        behreq(1500.0,	2.26,	48.89),
        behreq(1687.5,	3.03,	52.00),
        behreq(1875.0,	3.76,	54.74),
        behreq(2062.5,	4.43,	57.06),
        behreq(2250.0,	5.11,	58.45),
        behreq(2437.5,	5.71,	59.74),
        behreq(2625.0,	6.27,	60.62),
        behreq(2812.5,	6.79,	61.05),
        behreq(3000.0,	7.28,	61.49),
        behreq(3187.5,	7.75,	61.90),
        behreq(3375.0,	8.21,	62.24),
        behreq(3562.5,	8.62,	62.50),
        behreq(3750.0,	9.00,	62.72),
        behreq(3937.5,	9.37,	62.93),
        behreq(4125.0,	9.75,	63.04),
        behreq(4312.5,	10.08,	63.07),
        behreq(4500.0,	10.42,	63.04),
        behreq(4687.5,	10.76,	62.95),
        behreq(4875.0,	11.10,	62.74),
        behreq(5062.5,	11.39,	62.56),
        behreq(5250.0,	11.67,	62.38),
        behreq(5437.5,	11.95,	62.16),
        behreq(5625.0,	12.25,	62.10),
        behreq(5812.5,	12.48,	61.81),
        behreq(6000.0,	12.72,	61.62),
        behreq(6187.5,	12.95,	61.27),
        behreq(6375.0,	13.19,	61.14),
        behreq(6562.5,	13.35,	60.63),
        behreq(6750.0,	13.60,	60.27),
        behreq(6937.5,	13.79,	59.77),
        behreq(7125.0,	13.98,	59.45),
        behreq(7312.5,	14.16,	58.94),
        behreq(7500.0,	14.34,	58.67),
        behreq(7687.5,	14.49,	58.52),
        behreq(7875.0,	14.69,	58.29),
        behreq(8062.5,	14.85,	58.05),
        behreq(8250.0,	15.02,	57.58),
        behreq(8437.5,	15.18,	57.11),
        behreq(8625.0,	15.34,	56.70),
        behreq(8812.5,	15.45,	56.33),
        behreq(9000.0,	15.59,	55.98),
        behreq(9187.5,	15.72,	55.75),
        behreq(9375.0,	15.84,	55.45),
        behreq(9562.5,	15.96,	55.08),
        behreq(9750.0,	16.10,	54.68),
        behreq(9937.5,	16.23,	54.21),
        behreq(10125.0,	16.35,	53.87),
        behreq(10312.5,	16.47,	53.49),
        behreq(10500.0,	16.58,	53.15),
        behreq(10687.5,	16.69,	52.83),
        behreq(10875.0,	16.78,	52.42),
        behreq(11062.5,	16.90,	51.91),
        behreq(11250.0,	17.02,	51.64),
        behreq(11437.5,	17.11,	51.31),
        behreq(11625.0,	17.21,	50.94),
        behreq(11812.5,	17.32,	50.61),
        behreq(12000.0,	17.39,	50.33),
        behreq(12187.5,	17.47,	49.79),
        behreq(12375.0,	17.56,	49.48),
        behreq(12562.5,	17.63,	49.10),
        behreq(12750.0,	17.70,	48.79),
        behreq(12937.5,	17.77,	48.53),
        behreq(13125.0,	17.87,	48.37),
        behreq(13312.5,	17.97,	48.22),
        behreq(13500.0,	18.06,	47.86),
        behreq(13687.5,	18.15,	47.62),
        behreq(13875.0,	18.23,	47.23),
        behreq(14062.5,	18.27,	47.24),
        behreq(14250.0,	18.30,	46.64),
        behreq(14437.5,	18.39,	46.39),
        behreq(14625.0,	18.44,	46.04),
        behreq(14812.5,	18.52,	45.85),
        behreq(15000.0,	18.60,	45.22),
        behreq(15187.5,	18.65,	44.79),
        behreq(15375.0,	18.69,	44.47),
        behreq(15562.5,	18.76,	44.13),
        behreq(15750.0,	18.81,	43.75),
        behreq(15937.5,	18.87,	43.28),
        behreq(16125.0,	18.94,	43.15),
        behreq(16312.5,	18.99,	42.81),
        behreq(16500.0,	19.06,	42.55),
        behreq(16687.5,	19.07,	42.34),
        behreq(16875.0,	19.13,	42.14),
        behreq(17062.5,	19.19,	41.88),
        behreq(17250.0,	19.27,	41.68),
        behreq(17437.5,	19.31,	41.39),
        behreq(17625.0,	19.39,	41.15),
        behreq(17812.5,	19.43,	40.85),
        behreq(18000.0,	19.46,	40.53),
        behreq(18187.5,	19.45,	40.27),
        behreq(18375.0,	19.49,	39.96),
        behreq(18562.5,	19.53,	39.68),
        behreq(18750.0,	19.54,	39.43),
        behreq(18937.5,	19.57,	38.95),
        behreq(19125.0,	19.62,	38.66),
        behreq(19312.5,	19.66,	38.55),
        behreq(19500.0,	19.69,	38.35),
        behreq(19687.5,	19.74,	38.12),
        behreq(19875.0,	19.78,	38.06),
        behreq(20062.5,	19.83,	37.85),
        behreq(20250.0,	19.86,	37.45),
        behreq(20437.5,	19.92,	37.09),
        behreq(20625.0,	19.96,	37.03),
        behreq(20812.5,	20.00,	36.96),
        behreq(21000.0,	20.01,	36.64),
        behreq(21187.5,	20.06,	36.39),
        behreq(21375.0,	20.10,	36.64),
        behreq(21562.5,	20.13,	35.57),
        behreq(21750.0,	20.13,	25.88),
        behreq(21937.5,	20.16,	12.05),
        behreq(22125.0,	20.14,	-3.32),
        behreq(22312.5,	20.10,	-18.40),
        behreq(22500.0,	20.02,	-30.75),
        behreq(22687.5,	19.97,	-34.93),
        behreq(22875.0,	19.93,	-35.14),
        behreq(23062.5,	19.87,	-35.36),
        behreq(23250.0,	19.82,	-35.84),
        behreq(23437.5,	19.76,	-36.25),
        behreq(23625.0,	19.66,	-46.88),
        behreq(23812.5,	19.60,	-51.95),
        behreq(24000.0,	19.53,	-63.69) }};

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
// Frequency (Hz)	Magnitude (dB)	Phase (degrees)
#define gmleq(freq, ampdb, phdeg)  (ampdb)
t_EQStruct gMassenburgEQ = {
    48000, 256,
    {
        gmleq(0.0,	3.35,	0.00),
        gmleq(187.5,	3.32,	0.00),
        gmleq(375.0,	3.89,	0.00),
        gmleq(562.5,	4.60,	0.00),
        gmleq(750.0,	5.31,	0.00),
        gmleq(937.5,	6.34,	0.00),
        gmleq(1125.0,	7.09,	0.00),
        gmleq(1312.5,	7.88,	0.00),
        gmleq(1500.0,	8.62,	0.00),
        gmleq(1687.5,	9.30,	0.00),
        gmleq(1875.0,	9.96,	0.00),
        gmleq(2062.5,	10.36,	0.00),
        gmleq(2250.0,	10.92,	0.00),
        gmleq(2437.5,	11.36,	0.00),
        gmleq(2625.0,	11.78,	0.00),
        gmleq(2812.5,	12.09,	0.00),
        gmleq(3000.0,	12.48,	0.00),
        gmleq(3187.5,	12.76,	0.00),
        gmleq(3375.0,	13.04,	0.00),
        gmleq(3562.5,	13.27,	0.00),
        gmleq(3750.0,	13.49,	0.00),
        gmleq(3937.5,	13.70,	0.00),
        gmleq(4125.0,	13.92,	0.00),
        gmleq(4312.5,	14.08,	0.00),
        gmleq(4500.0,	14.24,	0.00),
        gmleq(4687.5,	14.39,	0.00),
        gmleq(4875.0,	14.50,	0.00),
        gmleq(5062.5,	14.64,	0.00),
        gmleq(5250.0,	14.75,	0.00),
        gmleq(5437.5,	14.85,	0.00),
        gmleq(5625.0,	14.95,	0.00),
        gmleq(5812.5,	15.03,	0.00),
        gmleq(6000.0,	15.12,	0.00),
        gmleq(6187.5,	15.22,	0.00),
        gmleq(6375.0,	15.29,	0.00),
        gmleq(6562.5,	15.33,	0.00),
        gmleq(6750.0,	15.41,	0.00),
        gmleq(6937.5,	15.46,	0.00),
        gmleq(7125.0,	15.51,	0.00),
        gmleq(7312.5,	15.57,	0.00),
        gmleq(7500.0,	15.63,	0.00),
        gmleq(7687.5,	15.67,	0.00),
        gmleq(7875.0,	15.70,	0.00),
        gmleq(8062.5,	15.74,	0.00),
        gmleq(8250.0,	15.80,	0.00),
        gmleq(8437.5,	15.82,	0.00),
        gmleq(8625.0,	15.84,	0.00),
        gmleq(8812.5,	15.88,	0.00),
        gmleq(9000.0,	15.92,	0.00),
        gmleq(9187.5,	15.93,	0.00),
        gmleq(9375.0,	15.95,	0.00),
        gmleq(9562.5,	15.99,	0.00),
        gmleq(9750.0,	16.01,	0.00),
        gmleq(9937.5,	16.01,	0.00),
        gmleq(10125.0,	16.03,	0.00),
        gmleq(10312.5,	16.06,	0.00),
        gmleq(10500.0,	16.08,	0.00),
        gmleq(10687.5,	16.08,	0.00),
        gmleq(10875.0,	16.11,	0.00),
        gmleq(11062.5,	16.13,	0.00),
        gmleq(11250.0,	16.13,	0.00),
        gmleq(11437.5,	16.14,	0.00),
        gmleq(11625.0,	16.16,	0.00),
        gmleq(11812.5,	16.18,	0.00),
        gmleq(12000.0,	16.17,	0.00),
        gmleq(12187.5,	16.18,	0.00),
        gmleq(12375.0,	16.20,	0.00),
        gmleq(12562.5,	16.20,	0.00),
        gmleq(12750.0,	16.19,	0.00),
        gmleq(12937.5,	16.21,	0.00),
        gmleq(13125.0,	16.23,	0.00),
        gmleq(13312.5,	16.21,	0.00),
        gmleq(13500.0,	16.21,	0.00),
        gmleq(13687.5,	16.23,	0.00),
        gmleq(13875.0,	16.24,	0.00),
        gmleq(14062.5,	16.23,	0.00),
        gmleq(14250.0,	16.23,	0.00),
        gmleq(14437.5,	16.25,	0.00),
        gmleq(14625.0,	16.25,	0.00),
        gmleq(14812.5,	16.23,	0.00),
        gmleq(15000.0,	16.26,	0.00),
        gmleq(15187.5,	16.28,	0.00),
        gmleq(15375.0,	16.26,	0.00),
        gmleq(15562.5,	16.24,	0.00),
        gmleq(15750.0,	16.27,	0.00),
        gmleq(15937.5,	16.27,	0.00),
        gmleq(16125.0,	16.24,	0.00),
        gmleq(16312.5,	16.26,	0.00),
        gmleq(16500.0,	16.27,	0.00),
        gmleq(16687.5,	16.27,	0.00),
        gmleq(16875.0,	16.25,	0.00),
        gmleq(17062.5,	16.26,	0.00),
        gmleq(17250.0,	16.28,	0.00),
        gmleq(17437.5,	16.25,	0.00),
        gmleq(17625.0,	16.25,	0.00),
        gmleq(17812.5,	16.27,	0.00),
        gmleq(18000.0,	16.27,	0.00),
        gmleq(18187.5,	16.24,	0.00),
        gmleq(18375.0,	16.24,	0.00),
        gmleq(18562.5,	16.25,	0.00),
        gmleq(18750.0,	16.26,	0.00),
        gmleq(18937.5,	16.24,	0.00),
        gmleq(19125.0,	16.26,	0.00),
        gmleq(19312.5,	16.26,	0.00),
        gmleq(19500.0,	16.25,	0.00),
        gmleq(19687.5,	16.24,	0.00),
        gmleq(19875.0,	16.28,	0.00),
        gmleq(20062.5,	16.27,	0.00),
        gmleq(20250.0,	16.24,	0.00),
        gmleq(20437.5,	16.25,	0.00),
        gmleq(20625.0,	16.28,	0.00),
        gmleq(20812.5,	16.26,	0.00),
        gmleq(21000.0,	16.24,	0.00),
        gmleq(21187.5,	16.27,	0.00),
        gmleq(21375.0,	16.29,	0.00),
        gmleq(21562.5,	16.26,	0.00),
        gmleq(21750.0,	16.26,	0.00),
        gmleq(21937.5,	16.32,	0.00),
        gmleq(22125.0,	16.36,	0.00),
        gmleq(22312.5,	16.45,	0.00),
        gmleq(22500.0,	16.75,	0.00),
        gmleq(22687.5,	17.13,	0.00),
        gmleq(22875.0,	17.69,	0.00),
        gmleq(23062.5,	18.56,	0.00),
        gmleq(23250.0,	19.73,	0.00),
        gmleq(23437.5,	21.17,	0.00),
        gmleq(23625.0,	23.42,	0.00),
        gmleq(23812.5,	24.49,	0.00),
        gmleq(24000.0,	23.96,	0.00) }};

// -------------------------------------------------------------
//

Float64 identity_Float64(Float64 v)
{
    return v;
}

void TCrescendo::set_postEQ(t_EQStruct *eqtbl, bool force)
{
    static Float64 postEQ[129];
    
    if((eqtbl != m_PostEQ_basis) || force)
    {
        interpolate_eqStruct(eqtbl, postEQ, &identity_Float64);
        m_PostEQ_basis = eqtbl;
        m_PostEQ = postEQ;
        invalidate_unified_filter();
    }
}

// -------------------------------------------------------------------------------------
void TCrescendo::set_postEQ(UInt32 ix)
{
	static t_EQStruct *eqps[] = {
		&gNullEQ,
		&gBehringerEQ,
		&gAH3300EQ,
		&gMassenburgEQ };
    
    if(ix > 3)
        ix = 3;
	set_postEQ(eqps[ix]);
}

// -------------------------------------------------------------------------------------

// Frequency (Hz)	Magnitude (dB)
#define DBM 0
#if DBM
#define hypercorr(freq, ampdb)  (ampdb)
#else
#define hypercorr(freq, ampdb)  (0.0)
#endif
// -2 dB peak dip @ 1.6 kHz Q = 0.707
// +3 dB hishelf @ 4 kHz Q = 0.707
// Converted from Phon space to dBSPL space
t_EQStruct gDMHyperCorrEQ = {
    48000, 256,
    {
        hypercorr(0.000,  -0.000),
        hypercorr(0.188,  -0.050),
        hypercorr(0.375,  -0.210),
        hypercorr(0.562,  -0.475),
        hypercorr(0.750,  -0.824),
        hypercorr(0.938,  -1.212),
        hypercorr(1.125,  -1.572),
        hypercorr(1.312,  -1.834),
        hypercorr(1.500,  -1.955),
        hypercorr(1.688,  -1.938),
        hypercorr(1.875,  -1.815),
        hypercorr(2.062,  -1.624),
        hypercorr(2.250,  -1.394),
        hypercorr(2.438,  -1.144),
        hypercorr(2.625,  -0.884),
        hypercorr(2.812,  -0.619),
        hypercorr(3.000,  -0.353),
        hypercorr(3.188,  -0.089),
        hypercorr(3.375,  0.170),
        hypercorr(3.562,  0.419),
        hypercorr(3.750,  0.656),
        hypercorr(3.938,  0.878),
        hypercorr(4.125,  1.082),
        hypercorr(4.312,  1.269),
        hypercorr(4.500,  1.439),
        hypercorr(4.688,  1.593),
        hypercorr(4.875,  1.732),
        hypercorr(5.062,  1.857),
        hypercorr(5.250,  1.969),
        hypercorr(5.438,  2.070),
        hypercorr(5.625,  2.161),
        hypercorr(5.812,  2.242),
        hypercorr(6.000,  2.314),
        hypercorr(6.188,  2.378),
        hypercorr(6.375,  2.434),
        hypercorr(6.562,  2.484),
        hypercorr(6.750,  2.528),
        hypercorr(6.938,  2.566),
        hypercorr(7.125,  2.599),
        hypercorr(7.312,  2.627),
        hypercorr(7.500,  2.652),
        hypercorr(7.688,  2.672),
        hypercorr(7.875,  2.689),
        hypercorr(8.062,  2.703),
        hypercorr(8.250,  2.714),
        hypercorr(8.438,  2.723),
        hypercorr(8.625,  2.728),
        hypercorr(8.812,  2.732),
        hypercorr(9.000,  2.733),
        hypercorr(9.188,  2.732),
        hypercorr(9.375,  2.730),
        hypercorr(9.562,  2.725),
        hypercorr(9.750,  2.718),
        hypercorr(9.938,  2.710),
        hypercorr(10.125,  2.699),
        hypercorr(10.312,  2.687),
        hypercorr(10.500,  2.674),
        hypercorr(10.688,  2.658),
        hypercorr(10.875,  2.641),
        hypercorr(11.062,  2.622),
        hypercorr(11.250,  2.602),
        hypercorr(11.438,  2.579),
        hypercorr(11.625,  2.555),
        hypercorr(11.812,  2.530),
        hypercorr(12.000,  2.502),
        hypercorr(12.188,  2.473),
        hypercorr(12.375,  2.442),
        hypercorr(12.562,  2.409),
        hypercorr(12.750,  2.374),
        hypercorr(12.938,  2.337),
        hypercorr(13.125,  2.298),
        hypercorr(13.312,  2.257),
        hypercorr(13.500,  2.215),
        hypercorr(13.688,  2.170),
        hypercorr(13.875,  2.123),
        hypercorr(14.062,  2.074),
        hypercorr(14.250,  2.023),
        hypercorr(14.438,  1.969),
        hypercorr(14.625,  1.913),
        hypercorr(14.812,  1.855),
        hypercorr(15.000,  1.794),
        hypercorr(15.188,  1.731),
        hypercorr(15.375,  1.666),
        hypercorr(15.562,  1.597),
        hypercorr(15.750,  1.527),
        hypercorr(15.938,  1.453),
        hypercorr(16.125,  1.377),
        hypercorr(16.312,  1.298),
        hypercorr(16.500,  1.216),
        hypercorr(16.688,  1.131),
        hypercorr(16.875,  1.043),
        hypercorr(17.062,  0.952),
        hypercorr(17.250,  0.858),
        hypercorr(17.438,  0.761),
        hypercorr(17.625,  0.661),
        hypercorr(17.812,  0.557),
        hypercorr(18.000,  0.450),
        hypercorr(18.188,  0.339),
        hypercorr(18.375,  0.225),
        hypercorr(18.562,  0.107),
        hypercorr(18.750,  0.000),
        hypercorr(18.938,  0.000),
        hypercorr(19.125,  0.000),
        hypercorr(19.312,  0.000),
        hypercorr(19.500,  0.000),
        hypercorr(19.688,  0.000),
        hypercorr(19.875,  0.000),
        hypercorr(20.062,  0.000),
        hypercorr(20.250,  0.000),
        hypercorr(20.438,  0.000),
        hypercorr(20.625,  0.000),
        hypercorr(20.812,  0.000),
        hypercorr(21.000,  0.000),
        hypercorr(21.188,  0.000),
        hypercorr(21.375,  0.000),
        hypercorr(21.562,  0.000),
        hypercorr(21.750,  0.000),
        hypercorr(21.938,  0.000),
        hypercorr(22.125,  0.000),
        hypercorr(22.312,  0.000),
        hypercorr(22.500,  0.000),
        hypercorr(22.688,  0.000),
        hypercorr(22.875,  0.000),
        hypercorr(23.062,  0.000),
        hypercorr(23.250,  0.000),
        hypercorr(23.438,  0.000),
        hypercorr(23.625,  0.000),
        hypercorr(23.812,  0.000),
        hypercorr(24.000,  0.000) }};


void TCrescendo::set_preEQ()
{
    static Float64 preEQ[129];
    static Float64 preEQAmpl[129];
    
    interpolate_eqStruct(&gDMHyperCorrEQ, preEQ, &identity_Float64, false);
    m_PreEQ = preEQ;
    interpolate_eqStruct(&gDMHyperCorrEQ, preEQAmpl, &ampl10, false);
    m_PreEQAmpl = preEQAmpl;
    invalidate_unified_filter();
}

