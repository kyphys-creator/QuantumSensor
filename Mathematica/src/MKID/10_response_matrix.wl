(* ::Package:: *)

(* ========================================================================== *)
(*  MKID - Response matrix (v_min integration of saved response functions)    *)
(*  TiN counterpart of TES/10_response_matrix.wl.                              *)
(*                                                                            *)
(*  Builds a response matrix from the response functions that 08 saved as      *)
(*  .wdx. The user gives a v_min range and a number of intervals; v_min is     *)
(*  split into that many equal intervals and each matrix element is the        *)
(*  integral of that bin's response over the interval:                        *)
(*                                                                            *)
(*       M[bin_i, interval_j] = Integral_{v_j}^{v_{j+1}} R_{bin_i}(v) dv       *)
(*                                                                            *)
(*  rows  = observed-energy bins,  columns = v_min intervals.                  *)
(*  Because the saved responses are order-1 (piecewise-linear) interpolating   *)
(*  functions, each integral is computed exactly by the trapezoid rule over    *)
(*  the function's own grid nodes inside the interval.                         *)
(*                                                                            *)
(*  Per-row effective window: each row (energy bin) is nonzero only above its   *)
(*  own kinematic threshold and then decays into a long ~6-10 decade tail. The  *)
(*  low edge (first nonzero) is kept exactly; the high tail is cut where the    *)
(*  row's cumulative integral reaches (1 - alpha) of its total (alpha default    *)
(*  0.001 -> keep 99.9% of every bin), zeroing entries beyond it. The matrix is  *)
(*  then trimmed to the populated column window: LEADING by the first row's      *)
(*  first nonzero (lowest threshold), TRAILING by the last row's last nonzero    *)
(*  (widest window); vmin.csv is trimmed at both ends to stay aligned.           *)
(*                                                                            *)
(*  This stage is light: it only imports the .wdx (self-contained) -- it does   *)
(*  NOT load the 01-06 pipeline.                                               *)
(*                                                                            *)
(*  Usage:                                                                     *)
(*      wolframscript -file 10_response_matrix.wl <name> <vminLo> <vminHi> <N> [alpha] *)
(*                                                                            *)
(*    <name>   base name of a .wdx in output/MKID/response_functions/, or ALL  *)
(*    <vminLo> <vminHi>  integration range [km/s] (clipped to the function's   *)
(*                       domain if it exceeds it)                              *)
(*    <N>      number of equal v_min intervals (matrix columns)               *)
(*    [alpha]  per-row tail-cut tolerance; keep (1-alpha) of each bin's        *)
(*             response (default 0.001). alpha = 0 disables the upper cut.      *)
(*                                                                            *)
(*  e.g.  wolframscript -file 10_response_matrix.wl TiN_q0M3_R5 1 800 1000      *)
(*        wolframscript -file 10_response_matrix.wl ALL 1 800 1000 0.0001      *)
(*                                                                            *)
(*  Output (in output/MKID/response_matrix/<mass>/<name>_v<lo>-<hi>_N<N>/):    *)
(*    matrix.csv   pure numeric matrix (nBins x N)                             *)
(*    vmin.csv     per-column {v_low, v_high, v_mid}                           *)
(*    bins.csv     per-row bin energy labels                                    *)
(* ========================================================================== *)

(* ============================ Arguments ============================ *)

commandLineArgs = Rest[$ScriptCommandLine];
If[Length[commandLineArgs] < 4,
  Print["usage: wolframscript -file 10_response_matrix.wl <name|ALL> <vminLo> <vminHi> <N> [alpha]"];
  Exit[1]];

nameArg    = commandLineArgs[[1]];
vminLoReq  = ToExpression[commandLineArgs[[2]]];
vminHiReq  = ToExpression[commandLineArgs[[3]]];
nIntervals = Round[ToExpression[commandLineArgs[[4]]]];
(* Per-row coverage tail-cut tolerance: keep (1 - alpha) of each bin's response
   (default 0.001 -> keep 99.9%). alpha = 0 disables the upper cut. *)
windowAlpha = If[Length[commandLineArgs] >= 5, ToExpression[commandLineArgs[[5]]], 0.001];
If[!(NumericQ[vminLoReq] && NumericQ[vminHiReq] && IntegerQ[nIntervals] && nIntervals >= 1
     && vminLoReq < vminHiReq && NumericQ[windowAlpha] && 0 <= windowAlpha < 1),
  Print["error: need numeric vminLo < vminHi, integer N >= 1, and 0 <= alpha < 1"];
  Exit[1]];


(* ============================ Setup ============================ *)

fileDir     = DirectoryName[$InputFileName];
mathDir     = ParentDirectory[ParentDirectory[fileDir]];
functionDir = FileNameJoin[{mathDir, "output", "MKID", "response_functions"}];
matrixDir   = FileNameJoin[{mathDir, "output", "MKID", "response_matrix"}];
If[!DirectoryQ[matrixDir], CreateDirectory[matrixDir, CreateIntermediateDirectories -> True]];

Get[FileNameJoin[{fileDir, "01_setup.wl"}]];


(* ============================ Integration ============================ *)

(* integrateIF[f, a, b] -> Integral_a^b f dv, exact for an order-1 (piecewise
   linear) InterpolatingFunction: trapezoid rule over f's grid nodes that lie
   inside (a, b), plus the endpoints a and b. *)
integrateIF[f_, a_, b_] := Module[{nodes, xs, ys},
  nodes = Flatten[f["Grid"]];
  xs = Join[{a}, Sort@Select[nodes, a < # < b &], {b}];
  ys = f /@ xs;
  Total[(Most[ys] + Rest[ys]) / 2 * Differences[xs]]
];


(* ============================ Build one matrix ============================ *)

extractMass[name_] := Module[{m},
  m = StringCases[name, "M" ~~ d : DigitCharacter .. :> "M" <> d];
  If[m === {}, "other", First[m]]
];

buildMatrix[wdxBaseName_] := Module[
  {wdxPath, data, responses, labels, fns, nBins, domLo, domHi,
   vminLo, vminHi, edges, mids, matrix, vminRows, firstCol, stem, massTag, outDir},

  wdxPath = FileNameJoin[{functionDir, wdxBaseName <> ".wdx"}];
  If[!FileExistsQ[wdxPath], Print["skip (not found): ", wdxBaseName <> ".wdx"]; Return[]];

  data      = Import[wdxPath];
  responses = data["responses"];
  labels    = Keys[responses];
  fns       = Values[responses];
  nBins     = Length[fns];

  (* clip the requested range to the interpolation domain *)
  {domLo, domHi} = First[fns]["Domain"][[1]];
  vminLo = Max[vminLoReq, domLo];
  vminHi = Min[vminHiReq, domHi];
  If[{vminLo, vminHi} =!= {vminLoReq, vminHiReq},
    Print["  note: range clipped to function domain {", domLo, ", ", domHi, "}"]];

  edges = Subdivide[vminLo, vminHi, nIntervals];          (* N+1 edges *)
  mids  = (Most[edges] + Rest[edges]) / 2;                 (* N midpoints *)

  (* M[i, j] = integral of bin-i response over interval j.
     The interpolation grid is in km/s; multiply by kps to convert dv to natural units. *)
  matrix   = kps Table[
    integrateIF[fns[[i]], edges[[j]], edges[[j + 1]]],
    {i, nBins}, {j, nIntervals}];
  vminRows = N /@ Transpose[{Most[edges], Rest[edges], mids}];

  (* ---- per-row effective window (zero the high-v tail) ----
     Each row rises sharply from a kinematic threshold (its first nonzero column
     -- the LOW edge is exact and kept as is) and then decays into a long
     ~6-10 decade tail. Keep the tail only up to where the row's cumulative
     integral reaches (1 - windowAlpha) of its total, zeroing everything beyond
     ("keep (1-alpha) of every bin's response"). windowAlpha = 0 disables it. *)
  rowHi[row_] := Module[{tot},
    tot = Total[row];
    If[tot <= 0, 0,
      Min[LengthWhile[Accumulate[row] / tot, # < 1 - windowAlpha &] + 1, Length[row]]]
  ];
  If[windowAlpha > 0,
    matrix = (PadRight[Take[#, rowHi[#]], nIntervals, 0.] & /@ matrix)];

  (* ---- trim to the column window actually populated ----
     LEADING: drop columns before the first nonzero of the FIRST row (lowest
     energy bin = lowest kinematic threshold = union left edge).
     TRAILING: drop columns after the last nonzero of the LAST row (highest
     energy bin = widest window = union right edge), scanned from the tail.
     v_min rows are trimmed identically so the columns stay aligned. *)
  firstCol = SelectFirst[Range[nIntervals], matrix[[1, #]] != 0 &, 1];
  lastCol  = SelectFirst[Range[nIntervals, 1, -1], matrix[[-1, #]] != 0 &, nIntervals];
  If[firstCol > 1 || lastCol < nIntervals,
    matrix   = matrix[[All, firstCol ;; lastCol]];
    vminRows = vminRows[[firstCol ;; lastCol]];
    Print["  window: kept cols ", firstCol, "..", lastCol, " / ", nIntervals,
          "  (", firstCol - 1, " leading + ", nIntervals - lastCol, " trailing trimmed)"]];

  stem    = wdxBaseName <> "_v" <> ToString[vminLoReq] <> "-" <> ToString[vminHiReq] <>
            "_N" <> ToString[nIntervals];
  massTag = extractMass[wdxBaseName];
  outDir  = FileNameJoin[{matrixDir, massTag, stem}];
  If[!DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];

  Export[FileNameJoin[{outDir, "matrix.csv"}], matrix];
  Export[FileNameJoin[{outDir, "vmin.csv"}], vminRows];
  Export[FileNameJoin[{outDir, "bins.csv"}],
    N /@ (ToExpression /@ StringCases[#, NumberString] & /@ labels)];

  Print["Saved: ", massTag, "/", stem, "/  (", nBins, " bins x ", Length[matrix[[1]]],
    " v_min intervals)"];
];


(* ============================ Dispatch ============================ *)

targets = If[ToUpperCase[nameArg] === "ALL",
  FileBaseName /@ FileNames["*.wdx", functionDir],
  {nameArg}];

If[targets === {},
  Print["no .wdx files in ", functionDir],
  buildMatrix /@ targets];
