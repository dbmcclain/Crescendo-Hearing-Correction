// crescendo_polys.h -- Polynomial approximations to the HC curves
//
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
// ------------------------------------------------------------------------
// Crescendo Nonlinear Compression Curves
// ------------------------------------------------------------------------------
// DBM/RAL 12/06
// weighted order (2,2) rational minimax fits
// {dbmin, p0,p1,p2,q1,q2}
//
// dbmin <= dbpwr <= 0
//
// gdb = (p0 + p1*x + p2*x^2)/(1.0 + q1*x + q2*x^2)
//

// ------------------------------------------------------------------------------
inline Float64 RANGE_REDUCE(Float64 dbpwr, Float64 *pcoffs)
{
  return (1.0 - 2.0*min(0.0, max(dbpwr-100.0, pcoffs[0]))/pcoffs[0]);
}

#if 0 // no need to limit gains on the low end
inline Float64 POST_SCALE(Float64 gdb)
{
  Float64 x = min(gdb, 24.0824);
  return ((x < 16.0) ? x : 16.0 + 0.75*(x-16.0));
}
#else
#define POST_SCALE(gdb)  (gdb)
#endif

inline Float64 APPROX(Float64 dbpwr, Float64 *pcoffs)
{
  Float64 x = RANGE_REDUCE(dbpwr, pcoffs);
  return ((pcoffs[3]*x+pcoffs[2])*x+pcoffs[1])/((pcoffs[5]*x+pcoffs[4])*x+1.0);
}

inline Float64 INTERPOLATE(Float64 dbpwr, Float64 frac, Float64 *pcoffs1, Float64 *pcoffs2)
{
  return ((0.0 == frac) 
	  ? APPROX(dbpwr, pcoffs1)
	  : (1.0-frac)*APPROX(dbpwr, pcoffs1) + frac*APPROX(dbpwr, pcoffs2));
}

// ------------------------------------------------------------------------------
#if GEN_POLY_TABLES
static Float64 gfit00[] = { -80.0,
			   0.0, 0.0, 0.0,
			   0.0, 0.0};
// (2,2) fits from NML
// unrestricted domain fits 20 dBSPL to 100 dBSPL
static Float64 gfit05[] = { -80.000000, 0.012153932845, -0.027972072658, 0.020234172436, 1.106218356485, 0.395151170592};
static Float64 gfit10[] = { -80.000000, 0.050195789141, -0.111617678367, 0.078211292363, 1.136573472741, 0.420641059795};
static Float64 gfit15[] = { -80.000000, 0.161335644299, -0.333041916426, 0.218357161984, 1.181995855308, 0.460606334581};
static Float64 gfit20[] = { -80.000000, 0.395090449522, -0.767315342693, 0.470643549457, 1.189576745240, 0.476954296066};
static Float64 gfit25[] = { -80.000000, 0.729531631616, -1.393947278250, 0.833516862035, 1.146557413417, 0.461619378482};
static Float64 gfit30[] = { -80.000000, 1.168306433997, -2.234778748792, 1.326413323836, 1.078196616528, 0.434322803072};
static Float64 gfit35[] = { -80.000000, 1.756814693160, -3.373932478990, 1.993065280096, 0.992069976176, 0.402011673342};
static Float64 gfit40[] = { -80.000000, 2.555796160098, -4.919994039434, 2.886270086186, 0.887762017845, 0.365618234421};
static Float64 gfit45[] = { -80.000000, 3.634032175130, -6.987899941065, 4.054920582166, 0.764540813057, 0.325257725544};
static Float64 gfit50[] = { -80.000000, 5.061765619606, -9.677313941142, 5.527783196897, 0.624093143413, 0.281860954703};
static Float64 gfit55[] = { -80.000000, 6.900453426629, -13.051866897679, 7.304307850380, 0.470569813788, 0.237132455000};
static Float64 gfit60[] = { -80.000000, 9.192293573511, -17.102019260069, 9.326713532491, 0.311587671121, 0.193715803621};
static Float64 gfit65[] = { -80.000000, 11.949868915007, -21.770617692341, 11.525134322195, 0.153794964781, 0.153476845029};
static Float64 gfit70[] = { -80.000000, 15.155610254241, -26.928113554838, 13.789208033938, 0.004699965880, 0.118201401328};
static Float64 gfit75[] = { -80.000000, 18.765714945524, -32.411640615103, 16.007719527784, -0.130068664218, 0.088779516867};
static Float64 gfit80[] = { -80.000000, 22.719979993526, -38.048453196509, 18.081404185623, -0.246995723140, 0.065333859211};
static Float64 gfit85[] = { -80.000000, 26.952568848659, -43.680698577395, 19.935687509170, -0.344735165018, 0.047387263449};
static Float64 gfit90[] = { -80.000000, 31.400741776793, -49.181772902765, 21.526279829902, -0.423735506702, 0.034112523180};
// ------------------------------------------------------------------------------

Float64* gfits[] = {
  gfit00, gfit05, gfit10, gfit15,
  gfit20, gfit25, gfit30, gfit35,
  gfit40, gfit45, gfit50, gfit55,
  gfit60, gfit65, gfit70, gfit75,
  gfit80, gfit85, gfit90, gfit00 };

#else
extern Float64* gfits[];
#endif // GEN_POLY_TABLES

// -- end of crescendo_polys.h -- //

