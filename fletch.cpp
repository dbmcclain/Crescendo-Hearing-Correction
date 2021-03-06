// fletch.cpp
// DM/RAL 10/07
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
#include "fletch.h"

// ----------------------------------------------------------
// Fletcher-Munson Conversion dBSPL <-> dBHL
//

// #define FMC(bark, fmdb)  (120.0/(120.0-(fmdb)))
#define FMC(bark, fmdb)  (fmdb)

#if 2 == NSUBBANDS // 50 Bark bands
Float64 gFletch[2*NFBANDS+1] = {
  FMC(0.0,36.607040586379284),
   FMC(0.5,19.583829825772057),
   FMC(1.0,13.218653604658158),
   FMC(1.5,9.801565253057348),
   FMC(2.0,7.640869063618908),
   FMC(2.5,6.13842123655003),
   FMC(3.0,5.026229587199509),
   FMC(3.5,4.165365617661266),
   FMC(4.0,3.476225862325933),
   FMC(4.5,2.8081075573208665),
   FMC(5.0,2.263677085859044),
   FMC(5.5,1.8086699922363607),
   FMC(6.0,1.3605851357049374),
   FMC(6.5,0.9781637061923902),
  FMC(7.0,0.6440728236675568),
  FMC(7.5,0.30551700965779105),
  FMC(8.0,0.0),
  FMC(8.5,-0.282878348047646),
  FMC(9.0,-0.5841109021291544),
  FMC(9.5,-0.9083913315623278),
  FMC(10.0,-1.231435106050327),
  FMC(10.5,-1.5950029204739167),
  FMC(11.0,-2.011036166528797),
  FMC(11.5,-2.4550644221857),
  FMC(12.0,-2.9731665146687902),
  FMC(12.5,-3.6204081077064387),
  FMC(13.0,-4.314276365717376),
  FMC(13.5,-5.136576775386715),
  FMC(14.0,-6.008521126685437),
  FMC(14.5,-6.908790891606429),
  FMC(15.0,-7.650313300924353),
  FMC(15.5,-8.22982331107552),
  FMC(16.0,-8.329078298628687),
  FMC(16.5,-7.808639602632553),
  FMC(17.0,-6.756611521473388),
  FMC(17.5,-5.026619609240782),
  FMC(18.0,-3.4854992195670462),
  FMC(18.5,-2.21098962422397),
  FMC(19.0,-1.498294321222133),
  FMC(19.5,-0.8872677584433148),
  FMC(20.0,-0.20242491073685187),
  FMC(20.5,0.8572431234816378),
  FMC(21.0,2.507996710880434),
  FMC(21.5,5.377062440067519),
  FMC(22.0,9.340813210043805),
  FMC(22.5,17.865538195724838),
  FMC(23.0,30.29976449959615),
  FMC(23.5,54.75728587300147),
  FMC(24.0,59.9),
  FMC(24.5,59.9),
  FMC(25.0,59.9)
};
#else  // 4 == NSUBBANDS //  100 Bark bands
Float64 gFletch[4*NFBANDS+1] = {
	FMC(0.0,44.87865326562556),
	FMC(0.25,31.775097624127902),
	FMC(0.5,24.54695253713337),
	FMC(0.75,19.90705961746405),
	FMC(1.0,16.649830058971624),
	FMC(1.25,14.223497469492678),
	FMC(1.5,12.33800978120197),
	FMC(1.75,10.825569649096503),
	FMC(2.0,9.582017027788663),
	FMC(2.25,8.539094896997448),
	FMC(2.5,7.650107225565556),
	FMC(2.75,6.881968859870822),
	FMC(3.0,6.210545157015508),
	FMC(3.25,5.617791302022448),
	FMC(3.5,5.089927137433246),
	FMC(3.75,4.616233681477599),
	FMC(4.0,4.188236684371947),
	FMC(4.25,3.79913881526647),
	FMC(4.5,3.443415997698662),
	FMC(4.75,3.116524755327591),
	FMC(5.0,2.8146862447166883),
	FMC(5.25,2.534724272770379),
	FMC(5.5,2.273941958599314),
	FMC(5.75,2.030026471202261),
	FMC(6.0,1.8009744310480796),
	FMC(6.25,1.5850326910697983),
	FMC(6.5,1.3806506706874076),
	FMC(6.75,1.186441431303193),
	FMC(7.0,1.0011493980848658),
	FMC(7.25,0.8236231453088156),
	FMC(7.5,0.6527920339342317),
	FMC(7.75,0.48764576299433715),
	FMC(8.0,0.32721610026447934),
	FMC(8.25,0.1705602133735491),
	FMC(8.5,0.01674514536674554),
	FMC(8.75,-0.1351669195548646),
	FMC(9.0,-0.2861328642148582),
	FMC(9.25,-0.4371424975707501),
	FMC(9.5,-0.5892321496611555),
	FMC(9.75,-0.7434979818010268),
	FMC(10.0,-0.9011088633281648),
	FMC(10.25,-1.0633185050117278),
	FMC(10.5,-1.2314763328681622),
	FMC(10.75,-1.407036320722393),
	FMC(11.0,-1.5915626598640596),
	FMC(11.25,-1.7867307129706858),
	FMC(11.5,-1.9943211608038356),
	FMC(11.75,-2.216204590663477),
	FMC(12.0,-2.4543129892592215),
	FMC(12.25,-2.7105936981998396),
	FMC(12.5,-2.986940402023723),
	FMC(12.75,-3.285094722603832),
	FMC(13.0,-3.6065111299360733),
	FMC(13.25,-3.9521773813499417),
	FMC(13.5,-4.322382932220708),
	FMC(13.75,-4.7164292534155425),
	FMC(14.0,-5.132279479295823),
	FMC(14.25,-5.56615124566545),
	FMC(14.5,-6.01206708450419),
	FMC(14.75,-6.46139247820069),
	FMC(15.0,-6.902413513348028),
	FMC(15.25,-7.320034026924512),
	FMC(15.5,-7.6957044398342305),
	FMC(15.75,-8.007726283662015),
	FMC(16.0,-8.232098294132932),
	FMC(16.25,-8.344066489588437),
	FMC(16.5,-8.32049048491474),
	FMC(16.75,-8.143016901817653),
	FMC(17.0,-7.8018387154723),
	FMC(17.25,-7.299518060308882),
	FMC(17.5,-6.65400417787484),
	FMC(17.75,-5.8996990442645405),
	FMC(18.0,-5.085394237217535),
	FMC(18.25,-4.268339360188142),
	FMC(18.5,-3.5047382841797164),
	FMC(18.75,-2.838473883331496),
	FMC(19.0,-2.2912800178701334),
	FMC(19.25,-1.8579886773264565),
	FMC(19.5,-1.5090615521658875),
	FMC(19.75,-1.1994680023367446),
	FMC(20.0,-0.8796497274713042),
	FMC(20.25,-0.5031587281375871),
	FMC(20.5,-0.027712748842096513),
	FMC(20.75,0.5897487466369511),
	FMC(21.0,1.4011539720550892),
	FMC(21.25,2.476556545401422),
	FMC(21.5,3.9151338348557494),
	FMC(21.75,5.86060489630621),
	FMC(22.0,8.525260506079322),
	FMC(22.25,12.229155625345241),
	FMC(22.5,17.465961035106336),
	FMC(22.75,25.017020320872664),
	FMC(23.0,36.155754858489246),
	FMC(23.25,53.02850717480218),
	FMC(23.5,79.39659829512134),
	FMC(23.75,80.0), // these were 120.0 but need to avoid 120.0
	FMC(24.0,80.0),  // to prevent division by zero
	FMC(24.25,80.0),
	FMC(24.5,80.0),
	FMC(24.75,80.0),
	FMC(25.0,80.0)
};
#endif // NSUBBANDS

// Performing ATH correction to dBHL by directly multiplying signal
// powers in the linear frequency FFT with an interpolated ATH filter
// helps to minimize seam artifacts in a swept tone as we cross from one
// analysis band to another.

// ATH in FFT domain at 48 kHz Fsamp
// FMC(freq Hz, mag dB)
t_EQStruct gATH = {
    48000, 256,
    {
        FMC(0.0,-120.0),
        FMC(187.5,10.501441764508407),
        FMC(375.0,4.570300728103248),
        FMC(562.5,2.326269900258007),
        FMC(750.0,1.081846054799787),
        FMC(937.5,0.23626622513140117),
        FMC(1125.0,-0.4351700862129402),
        FMC(1312.5,-1.0453238125790599),
        FMC(1500.0,-1.662699361301646),
        FMC(1687.5,-2.331709160732461),
        FMC(1875.0,-3.077416779401511),
        FMC(2062.5,-3.9045478590859717),
        FMC(2250.0,-4.79530126336542),
        FMC(2437.5,-5.708910865520167),
        FMC(2625.0,-6.584942764925055),
        FMC(2812.5,-7.351110723073354),
        FMC(3000.0,-7.934888344756869),
        FMC(3187.5,-8.276735814906893),
        FMC(3375.0,-8.341850723026962),
        FMC(3562.5,-8.127401022816783),
        FMC(3750.0,-7.6632780999731995),
        FMC(3937.5,-7.006185035895964),
        FMC(4125.0,-6.2287234042789486),
        FMC(4312.5,-5.4064063576973265),
        FMC(4500.0,-4.605800757446577),
        FMC(4687.5,-3.87628536341062),
        FMC(4875.0,-3.246571428989692),
        FMC(5062.5,-2.7257054110768486),
        FMC(5250.0,-2.307233025313461),
        FMC(5437.5,-1.9747861471083932),
        FMC(5625.0,-1.707524665516448),
        FMC(5812.5,-1.4844019704613045),
        FMC(6000.0,-1.2868500717836536),
        FMC(6187.5,-1.0999852228237664),
        FMC(6375.0,-0.9127194164264618),
        FMC(6562.5,-0.7172354825467648),
        FMC(6750.0,-0.508213393788052),
        FMC(6937.5,-0.2820650251323511),
        FMC(7125.0,-0.03630656799077103),
        FMC(7312.5,0.23089587431235925),
        FMC(7500.0,0.5210249340485178),
        FMC(7687.5,0.8353928144358371),
        FMC(7875.0,1.1752430665530853),
        FMC(8062.5,1.5418041159866083),
        FMC(8250.0,1.936315637485714),
        FMC(8437.5,2.3600407235907093),
        FMC(8625.0,2.814271107073623),
        FMC(8812.5,3.300329200932602),
        FMC(9000.0,3.819568772888758),
        FMC(9187.5,4.373375074824322),
        FMC(9375.0,4.963164774416102),
        FMC(9562.5,5.590385826553089),
        FMC(9750.0,6.2565173351275085),
        FMC(9937.5,6.963069421940318),
        FMC(10125.0,7.711583107154283),
        FMC(10312.5,8.503630201601867),
        FMC(10500.0,9.340813210043816),
        FMC(10687.5,10.224765244228016),
        FMC(10875.0,11.157149944641486),
        FMC(11062.5,12.139661409964086),
        FMC(11250.0,13.17402413335223),
        FMC(11437.5,14.261992944790844),
        FMC(11625.0,15.405352958847033),
        FMC(11812.5,16.605919527241447),
        FMC(12000.0,17.865538195724888),
        FMC(12187.5,19.18608466480861),
        FMC(12375.0,20.569464753951213),
        FMC(12562.5,22.017614368850328),
        FMC(12750.0,23.53249947152834),
        FMC(12937.5,25.116116052936278),
        FMC(13125.0,26.770490107831318),
        FMC(13312.5,28.49767761170891),
        FMC(13500.0,30.29976449959623),
        FMC(13687.5,32.17886664653279),
        FMC(13875.0,34.13712984958302),
        FMC(14062.5,36.176729811241614),
        FMC(14250.0,38.2998721241072),
        FMC(14437.5,40.50879225671206),
        FMC(14625.0,42.80575554040667),
        FMC(14812.5,45.19305715720841),
        FMC(15000.0,47.67302212853255),
        FMC(15187.5,50.248005304731066),
        FMC(15375.0,52.92039135537239),
        FMC(15562.5,55.69259476020069),
        FMC(15750.0,58.56705980072148),
        FMC(15937.5,61.546260552360174),
        FMC(16125.0,64.63270087715188),
        FMC(16312.5,67.82891441691814),
        FMC(16500.0,71.13746458689403),
        FMC(16687.5,74.56094456977131),
        FMC(16875.0,78.10197731012527),
        FMC(17062.5,81.76321550919855),
        FMC(17250.0,85.54734162001218),
        FMC(17437.5,89.45706784278366),
        FMC(17625.0,93.49513612062727),
        FMC(17812.5,97.6643181355168),
        FMC(18000.0,101.96741530449326),
        FMC(18187.5,106.40725877609921),
        FMC(18375.0,110.98670942702393),
        FMC(18562.5,115.70865785894628),
        FMC(18750.0,120.0),
        FMC(18937.5,120.0),
        FMC(19125.0,120.0),
        FMC(19312.5,120.0),
        FMC(19500.0,120.0),
        FMC(19687.5,120.0),
        FMC(19875.0,120.0),
        FMC(20062.5,120.0),
        FMC(20250.0,120.0),
        FMC(20437.5,120.0),
        FMC(20625.0,120.0),
        FMC(20812.5,120.0),
        FMC(21000.0,120.0),
        FMC(21187.5,120.0),
        FMC(21375.0,120.0),
        FMC(21562.5,120.0),
        FMC(21750.0,120.0),
        FMC(21937.5,120.0),
        FMC(22125.0,120.0),
        FMC(22312.5,120.0),
        FMC(22500.0,120.0),
        FMC(22687.5,120.0),
        FMC(22875.0,120.0),
        FMC(23062.5,120.0),
        FMC(23250.0,120.0),
        FMC(23437.5,120.0),
        FMC(23625.0,120.0),
        FMC(23812.5,120.0),
        FMC(24000.0,120.0)
    }
};

Float64 invAmpl10(Float64 v)
{
    return ampl10(-v);
}

void TCrescendo::compute_inverse_ATH_filter()
{
    interpolate_eqStruct(&gATH, m_InvATH, &invAmpl10);
    invalidate_unified_filter();
}


