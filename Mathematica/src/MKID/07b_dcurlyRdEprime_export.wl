(* ::Package:: *)

(* ========================================================================== *)
(*  MKID - d(curly R)/dE' kernel CSV export                                     *)
(*                                                                            *)
(*  Companion to 07_dcurlyRdEprime_plot.wl: instead of PDFs, write a plain CSV *)
(*  per mediator so the Python side can plot the differential response kernel  *)
(*  d curly R / d omega' (v_min) over the full v_min range.                     *)
(*                                                                            *)
(*  For each mediator (heavy n=0 / light n=2) it samples the summed left+right *)
(*  kinematic branches  KerRAll + KerRAlr  on a log-spaced v_min grid, for the  *)
(*  three DM masses (10/100/1000 MeV) at two observed energies E' (the bold    *)
(*  EBOLD and the faded 1 eV reference), in physical kg^-1 eV^-1.               *)
(*                                                                            *)
(*  Output columns:                                                            *)
(*    vmin, 10MeV_<EBOLD>eV, 10MeV_1eV, 100MeV_..., 1000MeV_...                 *)
(*                                                                            *)
(*  Usage:                                                                     *)
(*      wolframscript -file 07b_dcurlyRdEprime_export.wl  [nSamples]           *)
(*                                                                            *)
(*  Input:  06_response_defs.wl (kernel definitions)                           *)
(*  Output: output/MKID/dcurlyRdEprime_csv/TiN_{heavy,light}.csv                  *)
(* ========================================================================== *)

commandLineArgs = Rest[$ScriptCommandLine];
nSamples = If[Length[commandLineArgs] >= 1, Round[ToExpression[commandLineArgs[[1]]]], 600];

fileDir  = DirectoryName[$InputFileName];
mathDir  = ParentDirectory[ParentDirectory[fileDir]];
inputDir = FileNameJoin[{mathDir, "input", "MKID"}];
outDir   = FileNameJoin[{mathDir, "output", "MKID", "dcurlyRdEprime_csv"}];
If[!DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];

Get[FileNameJoin[{fileDir, "06_response_defs.wl"}]];
SetDirectory[inputDir];

(* ---- detector-specific knobs (see 07_dcurlyRdEprime_plot.wl) ---- *)
EBOLD   = 0.2;            (* MKID bold observed energy E' [eV] *)
sigPar  = MKIDsig;
kerL    = KerRTiNl;       (* left  kinematic branch *)
kerR    = KerRTiNr;       (* right kinematic branch *)
detTag  = "TiN";

toPhysicalKgEv = kg eV;  (* natural -> kg^-1 eV^-1 *)

(* summed kernel d curly R / d omega' at (mass md, mediator n, energy Ep[eV], v_min[km/s]) *)
kern[md_, n_, EpEv_, v_] := toPhysicalKgEv (
   kerL[md][n][EpEv eV, sigPar][v kps] + kerR[md][n][EpEv eV, sigPar][v kps]);

masses = {10, 100, 1000};                 (* MeV *)
energies = {EBOLD, 1};                     (* eV: bold, faded reference *)
vs = Exp[Subdivide[Log[1.], Log[10000.], nSamples - 1]];

header = Prepend[
  Flatten[Table[ToString[m] <> "MeV_" <> ToString[e] <> "eV", {m, masses}, {e, energies}]],
  "vmin"];

exportMediator[n_, label_] := Module[{table, outName},
  table = Table[
    Prepend[
      Flatten[Table[kern[m MeV, n, e, v], {m, masses}, {e, energies}]],
      v],
    {v, vs}];
  outName = detTag <> "_" <> label <> ".csv";
  Export[FileNameJoin[{outDir, outName}], Prepend[table, header]];
  Print["Saved: ", outName, "  (", Length[masses], " masses x ", Length[energies],
        " E', ", nSamples, " v_min pts)"];
];

Print["outDir = ", outDir];
exportMediator[0, "heavy"];
exportMediator[2, "light"];
