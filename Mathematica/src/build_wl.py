#!/usr/bin/env python3
"""Build TES/MKID .wl package files for 01-06 from scratch.

Plain Mathematica source — no Cell wrapping. Uses @@DET@@ as a placeholder
that is substituted with "TES" or "MKID" per detector, and @@HELPERS@@ for
the shared energySum/midpointSum/BoxPDF block in 06.
"""

from pathlib import Path

BASE = Path(__file__).parent


SETUP = r"""(* ========================================================================== *)
(*  @@DET@@ - Setup (Constants, Units, common Parameters)                       *)
(*  Natural units with all quantities expressed in GeV.                         *)
(* ========================================================================== *)

fileDir = DirectoryName[$InputFileName];
Direc = FileNameJoin[{ParentDirectory[ParentDirectory[fileDir]], "input", "@@DET@@"}];
SetDirectory[Direc];
Print["Direc = ", Direc];


(* ---- Energy base ---- *)
GeV = 10^9;
eV  = 10^(-9) GeV;
keV = 10^(-6) GeV;
MeV = 10^(-3) GeV;


(* ---- Other base units (GeV equivalents) ---- *)
grams = 5.62 * 10^23 GeV;
cm    = 1/(1.98 * 10^(-14) GeV);
sec   = 1/(6.58 * 10^(-25) GeV);
Kel   = 8.62 * 10^(-14) GeV;


(* ---- Mass scale ---- *)
ng = 10^(-9) grams;
\[Mu]g = 10^(-6) grams;
mg = 10^(-3) grams;
kg = 10^3 grams;
me = 0.5109989 MeV;
Mpl = 1.22 * 10^19 GeV;


(* ---- Length, Time ---- *)
km = 10^5 cm;
mtr = 10^2 cm;
mic = \[Mu]m = 10^(-4) cm;
nm = 10^(-7) cm;
yr = 365 * 24 * 3600 sec;
month = yr / 12;
day = 24 * 3600 sec;


(* ---- DM halo velocities (astro-ph/9710077v1, astro-ph/0611671) ---- *)
kps = km / sec;
v0 = 238 kps;
ve = 250 kps;
vesc = 544 kps;


(* ---- Particle physics constants and DM density ---- *)
alpha = 1 / 137;
\[Rho]DM = 0.4 GeV / cm^3;
\[Rho]b = 10^14 GeV / cm^3;
amukg = 1.660538782 * 10^(-27);


(* ---- Material densities (used by 04_material) ---- *)
rhoAl = 2.7 grams / cm^3;
rhoTiN = 5.4 grams / cm^3;
"""


FUNCTIONS_MATH = r"""(* ========================================================================== *)
(*  @@DET@@ - Math / Utility Functions                                          *)
(* ========================================================================== *)

fileDir = DirectoryName[$InputFileName];
Direc = FileNameJoin[{ParentDirectory[ParentDirectory[fileDir]], "input", "@@DET@@"}];
Get[FileNameJoin[{fileDir, "01_setup.wl"}]];
SetDirectory[Direc];
Print["Direc = ", Direc];


(* ============================ Math ============================ *)

(* LIntpl1D[x, y] — Build a piecewise-linear interpolation from parallel
   vectors x and y, returning an InterpolatingFunction (order 1). *)
LIntpl1D[x1_, y1_] := Module[{l1, Leng, lv},
  Leng = Length[x1];
  l1 = ConstantArray[0, Leng];
  For[lv = 1, lv <= Leng, lv++,
    l1[[lv]] = {{x1[[lv]]}, y1[[lv]]};
  ];
  ListInterpolation[l1, InterpolationOrder -> 1]
];

(* IntplArray[m] — Build an array of linear interpolations.
   m[[1]] is the common x grid; m[[k+1]] (k = 1..len-1) are y vectors. *)
IntplArray[m_] := Module[{len, intpl, lv},
  len = Length[m];
  intpl = ConstantArray[0, len - 1];
  For[lv = 1, lv <= len - 1, lv++,
    intpl[[lv]] = LIntpl1D[m[[1]], m[[lv + 1]]];
  ];
  intpl
];


(* ============================ Utility ============================ *)

(* Inword[word][size][x, y][color] — text annotation inset for figures.
   Places `word` at fractional position (x, y) of the image, with given
   font size and color, using a Times font. *)
Inword[word_][a_][x_, y_][color_] := Inset[
  Style[word, a, color, FontFamily -> "Times"],
  ImageScaled[{x, y}],
  {Center, Center}
];
"""


FUNCTIONS_RESPONSE = r"""(* ========================================================================== *)
(*  @@DET@@ - Domain / FDM / eta / Resolution functions                         *)
(* ========================================================================== *)

fileDir = DirectoryName[$InputFileName];
Direc = FileNameJoin[{ParentDirectory[ParentDirectory[fileDir]], "input", "@@DET@@"}];
Get[FileNameJoin[{fileDir, "02_functions_math.wl"}]];
SetDirectory[Direc];
Print["Direc = ", Direc];


(* ============================ Domain (kinematics) ============================ *)

(* Reduced mass of DM-target system *)
\[Mu]\[Chi]t[m\[Chi]_, mt_] := mt m\[Chi] / (mt + m\[Chi]);

(* Two kinematic branches of momentum transfer q at recoil energy Ee *)
ql[Ee_][m\[Chi]_][vmin_] := m\[Chi] vmin - Sqrt[m\[Chi]^2 vmin^2 - 2 m\[Chi] Ee];
qr[Ee_][m\[Chi]_][vmin_] := m\[Chi] vmin + Sqrt[m\[Chi]^2 vmin^2 - 2 m\[Chi] Ee];

(* Minimum q for given mass md (heavy mediator limit) *)
q0[Ee_][md_] := Sqrt[2 md Ee];

(* Jacobians dE/dq for left and right branches *)
Jacobl[Ee_][m\[Chi]_][vmin_] := -(m\[Chi] - m\[Chi]^2 vmin / Sqrt[m\[Chi]^2 vmin^2 - 2 m\[Chi] Ee]);
Jacobr[Ee_][m\[Chi]_][vmin_] :=  m\[Chi] + m\[Chi]^2 vmin / Sqrt[m\[Chi]^2 vmin^2 - 2 m\[Chi] Ee];

(* Minimum DM velocity *)
vmin[Ee_][m\[Chi]_][q_] := q / (2 m\[Chi]) + Ee / q;
vmin0[Ee_][md_] := Sqrt[2 Ee / md];

(* Energy-momentum relation for elastic scatter *)
Eq[q_][md_][v_] := q v - q^2 / (2 md);


(* ============================ FDM (DM form factor) ============================ *)

FDM[q_][n_] := (alpha me / q)^n;
FmedH[md_][q_] := ((md 240 kps)^2 + (100 md)^2) / (q^2 + (100 md)^2);
FmedL[md_][q_] := ((md 240 kps)^2 + 0) / (q^2 + 0);


(* ============================ eta (halo speed integrals) ============================ *)

KKf[v0_, vesc_] := v0^3 (-2. Exp[-vesc^2/v0^2] Pi (vesc/v0) + Pi^(3/2) Erf[vesc/v0]);

\[Eta]th[m\[Chi]_][vm_][v0_, ve_, vesc_] :=
  (\[Rho]DM \[Sigma]e / m\[Chi]) (v0^2 Pi / (2 ve KKf[v0, vesc])) Piecewise[{
    {-4 Exp[-vesc^2/v0^2] ve + Sqrt[Pi] v0 (Erf[(vm + ve)/v0] - Erf[(vm - ve)/v0]),
     vm < vesc - ve},
    {-2 Exp[-vesc^2/v0^2] (ve + vesc - vm) + Sqrt[Pi] v0 (Erf[vesc/v0] - Erf[(vm - ve)/v0]),
     vm < vesc + ve && vesc - ve < vm},
    {0, vm > vesc + ve}
  }];

\[Eta]td[m\[Chi]_][vm_][v0_, ve_, vesc_] :=
  (\[Rho]DM \[Sigma]e / m\[Chi]) (v0^2 Pi / (2 ve KKf[v0, vesc])) Piecewise[{
    {-4 Exp[-vesc^2/v0^2] ve + Sqrt[Pi] v0 (Erf[(vm + ve)/v0] - Erf[(vm - ve)/v0]),
     vm < vesc - ve},
    {-2 Exp[-vesc^2/v0^2] (ve + vesc - vm) + Sqrt[Pi] v0 (Erf[vesc/v0] - Erf[(vm - ve)/v0]),
     vm < vesc + ve && vesc - ve < vm},
    {0, vm > vesc + ve}
  }];


(* ============================ Resolution ============================ *)

sigE[Ep_, sig_] := sig Ep;
Eth[md_, vmax_] := md vmax^2 / 2;
"""


MATERIAL_TES = r"""(* ========================================================================== *)
(*  TES - Material (Al, Mermin model dielectric function)                       *)
(* ========================================================================== *)

fileDir = DirectoryName[$InputFileName];
Direc = FileNameJoin[{ParentDirectory[ParentDirectory[fileDir]], "input", "TES"}];
Get[FileNameJoin[{fileDir, "03_functions_response.wl"}]];
SetDirectory[Direc];
Print["Direc = ", Direc];


(* ============================ Al constants ============================ *)

rhoAl = 2.7 grams / cm^3;
mp = 0.94 GeV;
AAl = 27;
mAl = mp AAl;
a0Al = 4.046 * 10^(-8) cm;
vLA = 6.2 kps;


(* ============================ Al dielectric (Mermin model) ============================ *)

datAl = Import["Al_mermin.dat"];
datLAl = Length[datAl];

eps1Aldat  = ConstantArray[0, datLAl - 1];
eps2Aldat  = ConstantArray[0, datLAl - 1];
ImepsAldat = ConstantArray[0, datLAl - 1];

For[i = 2, i <= datLAl, i++,
  eps1Aldat[[i - 1]]  = {{datAl[[i, 1]], datAl[[i, 2]]}, datAl[[i, 3]]};
  eps2Aldat[[i - 1]]  = {{datAl[[i, 1]], datAl[[i, 2]]}, datAl[[i, 4]]};
  ImepsAldat[[i - 1]] = {{datAl[[i, 1]], datAl[[i, 2]]},
                         datAl[[i, 4]] / (datAl[[i, 4]]^2 + datAl[[i, 3]]^2)};
];

eps1AlIntpl  = Interpolation[eps1Aldat,  InterpolationOrder -> 1];
eps2AlIntpl  = Interpolation[eps2Aldat,  InterpolationOrder -> 1];
ImepsAlIntpl = Interpolation[ImepsAldat, InterpolationOrder -> 1];

ImepsAlf[w_, q_] := If[(0 <= w <= 99.3) && (0 <= q <= 37289.5), ImepsAlIntpl[w, q], 0];
eps1Alf[w_, q_]  := If[(0 <= w <= 99.3) && (0 <= q <= 37289.5), eps1AlIntpl[w, q],  0];
eps2Alf[w_, q_]  := If[(0 <= w <= 99.3) && (0 <= q <= 37289.5), eps2AlIntpl[w, q],  0];


(* ============================ Al phonon density of states ============================ *)

AlRawDos = Transpose[Delete[Import[FileNameJoin[{Direc, "Al_pDoS.dat"}]], 1]];
Print[{AlRawDos[[1, 1]], AlRawDos[[1, -1]]}];
DwAl = IntplArray[AlRawDos];
wDw = AlRawDos[[1, -1]];

Dw[w_?NumericQ] := Piecewise[{{DwAl[[1]][w], w <= wDw}, {0, w > wDw}}];

AveOmg    = NIntegrate[w     DwAl[[1]][w], {w, 0., wDw}];
AveInvOmg = NIntegrate[DwAl[[1]][w] / w,   {w, 0., wDw}];


(* ============================ Al form factor table ============================ *)

AlRawFn = Transpose[Delete[Import[FileNameJoin[{Direc, "Al_Fn.dat"}]], 1]];
Print[{AlRawFn[[1, 1]], AlRawFn[[1, -1]]}];
FnwAl = IntplArray[AlRawFn];

Fnw10GeV1[w_] := Piecewise[{{FnwAl[[1]][w], w <= 0.2}, {0, w > 0.2}}];
Fnw10GeV2[w_] := Piecewise[{{FnwAl[[2]][w], w <= 0.2}, {0, w > 0.2}}];
Fnw10GeV3[w_] := Piecewise[{{FnwAl[[3]][w], w <= 0.2}, {0, w > 0.2}}];
Fnw10GeV4[w_] := Piecewise[{{FnwAl[[4]][w], w <= 0.2}, {0, w > 0.2}}];
Fnw10GeV5[w_] := Piecewise[{{FnwAl[[5]][w], w <= 0.2}, {0, w > 0.2}}];
"""


MATERIAL_MKID = r"""(* ========================================================================== *)
(*  MKID - Material (TiN, analytic Lindhard dielectric function)                *)
(* ========================================================================== *)

fileDir = DirectoryName[$InputFileName];
Direc = FileNameJoin[{ParentDirectory[ParentDirectory[fileDir]], "input", "MKID"}];
Get[FileNameJoin[{fileDir, "03_functions_response.wl"}]];
SetDirectory[Direc];
Print["Direc = ", Direc];


(* ============================ TiN constants ============================ *)

rhoTiN = 5.4 grams / cm^3;
kFTiN = Sqrt[2 me (3.94 eV)];
vFTiN = kFTiN / me;
wpTiN = Sqrt[(4 Pi alpha / me) (2 me (3.94 eV))^(3/2) / (3 Pi^2)];
GammaTiN = 0.1 (3.94 eV);


(* ============================ TiN dielectric (Lindhard) ============================ *)

ImepsLTiN[w_, q_] := Im[-((1 + (3 wpTiN^2 / (q^2 vFTiN^2)) *
  (1/2 + (kFTiN / (4 q)) *
    (1 - ((q / (2 kFTiN)) - (w + I GammaTiN) / (q vFTiN))^2) *
    Log[((q / (2 kFTiN)) - (w + I GammaTiN) / (q vFTiN) + 1) /
        ((q / (2 kFTiN)) - (w + I GammaTiN) / (q vFTiN) - 1)] +
    (kFTiN / (4 q)) *
    (1 - ((q / (2 kFTiN)) + (w + I GammaTiN) / (q vFTiN))^2) *
    Log[((q / (2 kFTiN)) + (w + I GammaTiN) / (q vFTiN) + 1) /
        ((q / (2 kFTiN)) + (w + I GammaTiN) / (q vFTiN) - 1)]))^(-1))];
"""


PARAMETERS = r"""(* ========================================================================== *)
(*  @@DET@@ - Calculation Parameters (cross sections, resolutions, exposure)    *)
(* ========================================================================== *)

fileDir = DirectoryName[$InputFileName];
Direc = FileNameJoin[{ParentDirectory[ParentDirectory[fileDir]], "input", "@@DET@@"}];
Get[FileNameJoin[{fileDir, "04_material.wl"}]];
SetDirectory[Direc];
Print["Direc = ", Direc];


(* ---- Reference cross sections ---- *)
\[Sigma]e = 10^(-30) cm^2;
\[Sigma]N = 10^(-32) cm^2;


(* ---- Energy resolutions ---- *)
TESsig = (0.068 / 0.8) / 2.355;
MKIDsig = 0.3 / 2.355;


(* ---- Detector exposures (active mass times time) ---- *)
TiNexp = 10^7 rhoTiN (1 \[Mu]m) (50 \[Mu]m) (22 nm) month;
Alexp = (8200 ng) month;
"""


RESPONSE_HELPERS = r"""(* ============================ Helpers ============================ *)

(* energySum[f, E1, E2, dE]
   observed energy w Riemann sum: (sum_i f(w_i)) * dE * eV.
   - f is a 1-arg pure function
   - dE optional (default 0.01) *)
energySum[f_, E1_, E2_, dE_:0.01] := dE eV Total[f /@ Range[E1, E2, dE]];

(* midpointSum[f, vmin, vmax, n]
   velocity midpoint-rule integral over n midpoints in [vmin, vmax]. *)
midpointSum[f_, vmin_, vmax_, n_Integer] := With[
  {dv = (vmax - vmin)/n},
  dv Total[f /@ Table[vmin + (i - 0.5) dv, {i, n}]]
];


(* ---- Resolution ---- *)
(* BoxPDF[Ep, sig][ER] — rectangular PDF centred at Ep with half-width sig. *)
BoxPDF[Ep_, sig_][ER_] := UnitStep[sig - Abs[ER - Ep]] / (2 sig);
"""


RESPONSE_TES = r"""(* ========================================================================== *)
(*  TES - Response Function Definitions                                         *)
(* ========================================================================== *)

fileDir = DirectoryName[$InputFileName];
Direc = FileNameJoin[{ParentDirectory[ParentDirectory[fileDir]], "input", "TES"}];
Get[FileNameJoin[{fileDir, "05_parameters.wl"}]];
SetDirectory[Direc];
Print["Direc = ", Direc];


@@HELPERS@@


(* ============================ Raw differential rate kernels (Al) ============================ *)

kerRAll[E_][m\[Chi]_][n_][vmin_] := With[
  {q = ql[E][m\[Chi]][vmin]},
  If[Internal`RealValuedNumericQ[q],
    (1/rhoAl) (1/(8 \[Pi]^2 alpha)) (1/(\[Mu]\[Chi]t[m\[Chi], me])^2) *
      Jacobl[E][m\[Chi]][vmin] q^3 *
      FDM[q][n]^2 *
      ImepsAlf[E, q],
    0
  ]
];

kerRAlr[E_][m\[Chi]_][n_][vmin_] := With[
  {q = qr[E][m\[Chi]][vmin]},
  If[Internal`RealValuedNumericQ[q],
    (1/rhoAl) (1/(8 \[Pi]^2 alpha)) (1/(\[Mu]\[Chi]t[m\[Chi], me])^2) *
      Jacobr[E][m\[Chi]][vmin] q^3 *
      FDM[q][n]^2 *
      ImepsAlf[E, q],
    0
  ]
];


(* ============================ Integration domain ============================ *)

RangeAlGen[vmin_][md_][Ep_][Sig_] := Module[
  {emin, emax},
  emax = Min[(md vmin^2)/2, Ep + Sig];
  emin = Max[0.1 eV, Ep - Sig];
  If[emin < emax, {emin, emax}, {0, 0}]
];


(* ============================ Resolution-convolved kernels ============================ *)

IntkerRAll[E1_, E2_][m\[Chi]_][n_][Ep_, sig_][vmin_] :=
  energySum[
    PDF[NormalDistribution[Ep, sig], #] kerRAll[#][m\[Chi]][n][vmin] &,
    E1, E2, 0.001
  ];

IntkerRAlr[E1_, E2_][m\[Chi]_][n_][Ep_, sig_][vmin_] :=
  energySum[
    PDF[NormalDistribution[Ep, sig], #] kerRAlr[#][m\[Chi]][n][vmin] &,
    E1, E2, 0.001
  ];


(* ============================ Auto-range user kernels ============================ *)

KerRAll[md_][n_][Ep_, sig_][vmin_] := Module[
  {r1, r2},
  {r1, r2} = RangeAlGen[vmin][md][Ep][sig];
  If[r1 == 0 && r2 == 0, 0, IntkerRAll[r1, r2][md][n][Ep, sig][vmin]]
];

KerRAlr[md_][n_][Ep_, sig_][vmin_] := Module[
  {r1, r2},
  {r1, r2} = RangeAlGen[vmin][md][Ep][sig];
  If[r1 == 0 && r2 == 0, 0, IntkerRAlr[r1, r2][md][n][Ep, sig][vmin]]
];


(* ============================ Continuous-vm response ============================ *)

CRTESLeft[md_, n_][E1_, E2_][vm_?NumericQ] :=
  energySum[KerRAll[md][n][#, TESsig #][vm kps] &, E1, E2];

CRTESRight[md_, n_][E1_, E2_][vm_?NumericQ] :=
  energySum[KerRAlr[md][n][#, TESsig #][vm kps] &, E1, E2];

CRTES[md_, n_][E1_, E2_][vm_?NumericQ] :=
  energySum[(KerRAll[md][n][#, TESsig #][vm kps] +
             KerRAlr[md][n][#, TESsig #][vm kps]) &, E1, E2];


(* ============================ Velocity-integrated response ============================ *)

CRintTESLeft[md_, n_][E1_, E2_][vmin_?NumericQ, vmax_?NumericQ, ns_Integer] :=
  midpointSum[CRTESLeft[md, n][E1, E2][#] &, vmin, vmax, ns];

CRintTESRight[md_, n_][E1_, E2_][vmin_?NumericQ, vmax_?NumericQ, ns_Integer] :=
  midpointSum[CRTESRight[md, n][E1, E2][#] &, vmin, vmax, ns];

CRintTES[md_, n_][E1_, E2_][vmin_?NumericQ, vmax_?NumericQ, ns_Integer] :=
  midpointSum[CRTES[md, n][E1, E2][#] &, vmin, vmax, ns];

CRintTESWeighted[md_, n_][E1_, E2_][vmin_?NumericQ, vmax_?NumericQ, ns_Integer, eta_] :=
  midpointSum[CRTES[md, n][E1, E2][#] eta[#] &, vmin, vmax, ns];


(* ============================ Legacy table API ============================ *)

CRvmTES[md_, n_][E1_, E2_] :=
  Table[CRTES[md, n][E1, E2][vm], {vm, 5, 108, 0.51758}];
CRvmTESLeft[md_, n_][E1_, E2_] :=
  Table[CRTESLeft[md, n][E1, E2][vm], {vm, 5, 200, 196/800}];
CRvmTESRight[md_, n_][E1_, E2_] :=
  Table[CRTESRight[md, n][E1, E2][vm], {vm, 5, 200, 196/800}];
"""


RESPONSE_MKID = r"""(* ========================================================================== *)
(*  MKID - Response Function Definitions (TiN)                                  *)
(* ========================================================================== *)

fileDir = DirectoryName[$InputFileName];
Direc = FileNameJoin[{ParentDirectory[ParentDirectory[fileDir]], "input", "MKID"}];
Get[FileNameJoin[{fileDir, "05_parameters.wl"}]];
SetDirectory[Direc];
Print["Direc = ", Direc];


@@HELPERS@@


(* ============================ Raw differential rate kernels (TiN) ============================ *)

kerRTiNl[E_][m\[Chi]_][n_][vmin_] := With[
  {q = ql[E][m\[Chi]][vmin]},
  If[Internal`RealValuedNumericQ[q],
    (1/rhoTiN) (1/(8 \[Pi]^2 alpha)) (1/(\[Mu]\[Chi]t[m\[Chi], me])^2) *
      Jacobl[E][m\[Chi]][vmin] q^3 *
      FDM[q][n]^2 *
      ImepsLTiN[E, q],
    0
  ]
];

kerRTiNr[E_][m\[Chi]_][n_][vmin_] := With[
  {q = qr[E][m\[Chi]][vmin]},
  If[Internal`RealValuedNumericQ[q],
    (1/rhoTiN) (1/(8 \[Pi]^2 alpha)) (1/(\[Mu]\[Chi]t[m\[Chi], me])^2) *
      Jacobr[E][m\[Chi]][vmin] q^3 *
      FDM[q][n]^2 *
      ImepsLTiN[E, q],
    0
  ]
];


(* ============================ Integration domain ============================ *)

RangeTiNGen[vmin_][md_][Ep_][Sig_] := Module[
  {emin, emax},
  emax = Min[(md vmin^2)/2, Ep + Sig];
  emin = Max[0.1 eV, Ep - Sig];
  If[emin < emax, {emin, emax}, {0, 0}]
];


(* ============================ Resolution-convolved kernels ============================ *)

IntkerRTiNl[E1_, E2_][m\[Chi]_][n_][Ep_, sig_][vmin_] :=
  energySum[
    PDF[NormalDistribution[Ep, sig], #] kerRTiNl[#][m\[Chi]][n][vmin] &,
    E1, E2, 0.001
  ];

IntkerRTiNr[E1_, E2_][m\[Chi]_][n_][Ep_, sig_][vmin_] :=
  energySum[
    PDF[NormalDistribution[Ep, sig], #] kerRTiNr[#][m\[Chi]][n][vmin] &,
    E1, E2, 0.001
  ];


(* ============================ Auto-range user kernels ============================ *)

KerRTiNl[md_][n_][Ep_, sig_][vmin_] := Module[
  {r1, r2},
  {r1, r2} = RangeTiNGen[vmin][md][Ep][sig];
  If[r1 == 0 && r2 == 0, 0, IntkerRTiNl[r1, r2][md][n][Ep, sig][vmin]]
];

KerRTiNr[md_][n_][Ep_, sig_][vmin_] := Module[
  {r1, r2},
  {r1, r2} = RangeTiNGen[vmin][md][Ep][sig];
  If[r1 == 0 && r2 == 0, 0, IntkerRTiNr[r1, r2][md][n][Ep, sig][vmin]]
];


(* ============================ Continuous-vm response ============================ *)

CRTiNLeft[md_, n_][E1_, E2_][vm_?NumericQ] :=
  energySum[KerRTiNl[md][n][#, MKIDsig #][vm kps] &, E1, E2];

CRTiNRight[md_, n_][E1_, E2_][vm_?NumericQ] :=
  energySum[KerRTiNr[md][n][#, MKIDsig #][vm kps] &, E1, E2];

CRTiN[md_, n_][E1_, E2_][vm_?NumericQ] :=
  energySum[(KerRTiNl[md][n][#, MKIDsig #][vm kps] +
             KerRTiNr[md][n][#, MKIDsig #][vm kps]) &, E1, E2];


(* ============================ Velocity-integrated response ============================ *)

CRintTiNLeft[md_, n_][E1_, E2_][vmin_?NumericQ, vmax_?NumericQ, ns_Integer] :=
  midpointSum[CRTiNLeft[md, n][E1, E2][#] &, vmin, vmax, ns];

CRintTiNRight[md_, n_][E1_, E2_][vmin_?NumericQ, vmax_?NumericQ, ns_Integer] :=
  midpointSum[CRTiNRight[md, n][E1, E2][#] &, vmin, vmax, ns];

CRintTiN[md_, n_][E1_, E2_][vmin_?NumericQ, vmax_?NumericQ, ns_Integer] :=
  midpointSum[CRTiN[md, n][E1, E2][#] &, vmin, vmax, ns];

CRintTiNWeighted[md_, n_][E1_, E2_][vmin_?NumericQ, vmax_?NumericQ, ns_Integer, eta_] :=
  midpointSum[CRTiN[md, n][E1, E2][#] eta[#] &, vmin, vmax, ns];


(* ============================ Legacy table API ============================ *)

CRvmTiN[md_, n_][E1_, E2_] :=
  Table[CRTiN[md, n][E1, E2][vm], {vm, 6, 108, 0.51256}];
CRvmTiNLeft[md_, n_][E1_, E2_] :=
  Table[CRTiNLeft[md, n][E1, E2][vm], {vm, 6, 200, 195/800}];
CRvmTiNRight[md_, n_][E1_, E2_] :=
  Table[CRTiNRight[md, n][E1, E2][vm], {vm, 6, 200, 195/800}];
"""


def render(template: str, det: str) -> str:
    return template.replace("@@DET@@", det).replace("@@HELPERS@@", RESPONSE_HELPERS.strip())


def main():
    files = {
        "TES": {
            "01_setup.wl":              render(SETUP, "TES"),
            "02_functions_math.wl":     render(FUNCTIONS_MATH, "TES"),
            "03_functions_response.wl": render(FUNCTIONS_RESPONSE, "TES"),
            "04_material.wl":           MATERIAL_TES,
            "05_parameters.wl":         render(PARAMETERS, "TES"),
            "06_response_defs.wl":      render(RESPONSE_TES, "TES"),
        },
        "MKID": {
            "01_setup.wl":              render(SETUP, "MKID"),
            "02_functions_math.wl":     render(FUNCTIONS_MATH, "MKID"),
            "03_functions_response.wl": render(FUNCTIONS_RESPONSE, "MKID"),
            "04_material.wl":           MATERIAL_MKID,
            "05_parameters.wl":         render(PARAMETERS, "MKID"),
            "06_response_defs.wl":      render(RESPONSE_MKID, "MKID"),
        },
    }
    for det, group in files.items():
        outdir = BASE / det
        print(f"== {det} ==")
        for name, content in group.items():
            out = outdir / name
            out.write_text(content)
            print(f"  wrote {out.name}  ({len(content)} bytes)")


if __name__ == "__main__":
    main()
