(* ::Package:: *)

(* ========================================================================== *)
(*  TES - Response matrix (v_min integration of saved response functions)     *)
(*  (formerly 10_data; repurposed.)                                           *)
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
(*  This stage is light: it only imports the .wdx (self-contained) -- it does   *)
(*  NOT load the 01-06 pipeline.                                               *)
(*                                                                            *)
(*  Usage:                                                                     *)
(*      wolframscript -file 10_response_matrix.wl <name> <vminLo> <vminHi> <N> *)
(*                                                                            *)
(*    <name>   base name of a .wdx in output/TES/response_functions/, or ALL   *)
(*    <vminLo> <vminHi>  integration range [km/s] (clipped to the function's   *)
(*                       domain if it exceeds it)                              *)
(*    <N>      number of equal v_min intervals (matrix columns)               *)
(*                                                                            *)
(*  e.g.  wolframscript -file 10_response_matrix.wl Al_q0M3_R5 5 800 200       *)
(*        wolframscript -file 10_response_matrix.wl ALL 5 500 100             *)
(*                                                                            *)
(*  Output (in output/TES/response_matrix/):                                   *)
(*    <name>_v<lo>-<hi>_N<N>.csv        pure numeric matrix (nBins x N)        *)
(*    <name>_v<lo>-<hi>_N<N>_vmin.csv   per-column {v_low, v_high, v_mid}      *)
(*    <name>_v<lo>-<hi>_N<N>_bins.csv   per-row bin energy labels              *)
(* ========================================================================== *)

(* ============================ Arguments ============================ *)

commandLineArgs = Rest[$ScriptCommandLine];
If[Length[commandLineArgs] < 4,
  Print["usage: wolframscript -file 10_response_matrix.wl <name|ALL> <vminLo> <vminHi> <N>"];
  Exit[1]];

nameArg    = commandLineArgs[[1]];
vminLoReq  = ToExpression[commandLineArgs[[2]]];
vminHiReq  = ToExpression[commandLineArgs[[3]]];
nIntervals = Round[ToExpression[commandLineArgs[[4]]]];
If[!(NumericQ[vminLoReq] && NumericQ[vminHiReq] && IntegerQ[nIntervals] && nIntervals >= 1
     && vminLoReq < vminHiReq),
  Print["error: need numeric vminLo < vminHi and integer N >= 1"];
  Exit[1]];


(* ============================ Setup ============================ *)

fileDir     = DirectoryName[$InputFileName];
mathDir     = ParentDirectory[ParentDirectory[fileDir]];
functionDir = FileNameJoin[{mathDir, "output", "TES", "response_functions"}];
matrixDir   = FileNameJoin[{mathDir, "output", "TES", "response_matrix"}];
If[!DirectoryQ[matrixDir], CreateDirectory[matrixDir, CreateIntermediateDirectories -> True]];


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

buildMatrix[wdxBaseName_] := Module[
  {wdxPath, data, responses, labels, fns, nBins, domLo, domHi,
   vminLo, vminHi, edges, mids, matrix, stem},

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

  (* M[i, j] = integral of bin-i response over interval j *)
  matrix = Table[
    integrateIF[fns[[i]], edges[[j]], edges[[j + 1]]],
    {i, nBins}, {j, nIntervals}];

  stem = wdxBaseName <> "_v" <> ToString[vminLoReq] <> "-" <> ToString[vminHiReq] <>
         "_N" <> ToString[nIntervals];
  Export[FileNameJoin[{matrixDir, stem <> ".csv"}], matrix];
  Export[FileNameJoin[{matrixDir, stem <> "_vmin.csv"}],
    Transpose[{Most[edges], Rest[edges], mids}]];
  Export[FileNameJoin[{matrixDir, stem <> "_bins.csv"}], List /@ labels];

  Print["Saved: ", stem, ".csv  (", nBins, " bins x ", nIntervals, " v_min intervals)"];
];


(* ============================ Dispatch ============================ *)

targets = If[ToUpperCase[nameArg] === "ALL",
  FileBaseName /@ FileNames["*.wdx", functionDir],
  {nameArg}];

If[targets === {},
  Print["no .wdx files in ", functionDir],
  buildMatrix /@ targets];
