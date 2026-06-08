(* ========================================================================== *)
(*  MKID - Response Function Definitions (TiN)                                  *)
(* ========================================================================== *)

fileDir = DirectoryName[$InputFileName];
Direc = FileNameJoin[{ParentDirectory[ParentDirectory[fileDir]], "input", "MKID"}];
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
