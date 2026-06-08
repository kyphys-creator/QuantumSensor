(* ========================================================================== *)
(*  TES - Response Function Definitions                                         *)
(* ========================================================================== *)

fileDir = DirectoryName[$InputFileName];
Direc = FileNameJoin[{ParentDirectory[ParentDirectory[fileDir]], "input", "TES"}];
Get[FileNameJoin[{fileDir, "05_parameters.wl"}]];
SetDirectory[Direc];
Print["Direc = ", Direc];


(* ============================ Helpers ============================ *)

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
