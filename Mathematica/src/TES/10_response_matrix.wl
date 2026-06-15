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
(*  Per-row central window (two-sided, highest-density): each row peaks near    *)
(*  its kinematic edge and decays into a long tail. Keep only the highest-       *)
(*  density contiguous region around the peak covering a fraction `alpha` of     *)
(*  the row's area, zeroing BOTH the low-v_min rise and the high-v_min tail.     *)
(*  This localises the response so the monotone inverse recovers eta to the      *)
(*  window edge. The matrix is then trimmed to the union of the per-row windows; *)
(*  vmin.csv is trimmed identically so the columns stay aligned.                 *)
(*                                                                            *)
(*  This stage is light: it only imports the .wdx (self-contained) -- it does   *)
(*  NOT load the 01-06 pipeline.                                               *)
(*                                                                            *)
(*  Usage:                                                                     *)
(*      wolframscript -file 10_response_matrix.wl <name> <vminLo> <vminHi> <N> [alpha] *)
(*                                                                            *)
(*    <name>   base name of a .wdx in output/TES/response_functions/, or ALL   *)
(*    <vminLo> <vminHi>  integration range [km/s] (clipped to the function's   *)
(*                       domain if it exceeds it)                              *)
(*    <N>      number of equal v_min intervals (matrix columns)               *)
(*    [alpha]  central area fraction kept per row around its peak (two-sided;  *)
(*             default 0.5 = keep 50%). alpha = 0 disables the cut (full row).  *)
(*                                                                            *)
(*  e.g.  wolframscript -file 10_response_matrix.wl Al_q0M3_R5 1 800 1000       *)
(*        wolframscript -file 10_response_matrix.wl ALL 1 800 1000 0.0001      *)
(*                                                                            *)
(*  Output (in output/TES/response_matrix/<mass>/<name>_v<lo>-<hi>_N<N>/):     *)
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
(* Central two-sided area cut: keep the highest-density region around each
   bin's peak covering fraction `alpha` of its area, zeroing both the low-v_min
   rise and the high-v_min tail (default 0.5 -> keep 50%). Localising the
   response this way lets the monotone inverse recover eta to the window edge
   (validated in sandbox/vertex_weight_test). alpha = 0 disables the cut. *)
windowAlpha = If[Length[commandLineArgs] >= 5, ToExpression[commandLineArgs[[5]]], 0.5];
If[!(NumericQ[vminLoReq] && NumericQ[vminHiReq] && IntegerQ[nIntervals] && nIntervals >= 1
     && vminLoReq < vminHiReq && NumericQ[windowAlpha] && 0 <= windowAlpha < 1),
  Print["error: need numeric vminLo < vminHi, integer N >= 1, and 0 <= alpha < 1"];
  Exit[1]];


(* ============================ Setup ============================ *)

fileDir     = DirectoryName[$InputFileName];
mathDir     = ParentDirectory[ParentDirectory[fileDir]];
functionDir = FileNameJoin[{mathDir, "output", "TES", "response_functions"}];
matrixDir   = FileNameJoin[{mathDir, "output", "TES", "response_matrix"}];
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

  (* ---- per-row central window (two-sided, highest-density) ----
     Each matrix row M[i,j] is the response integrated over interval j, so every
     cell is an area element. Grow a window outward from the row's peak cell,
     each step annexing the larger neighbouring cell, until the kept area reaches
     windowAlpha of the row's total; zero everything outside (both the low-v_min
     rise and the high-v_min tail). windowAlpha = 0 keeps the full row. *)
  centralWindow[row_] := Module[{total, n, pk, lo, hi, acc, lv, rv},
    total = Total[row];
    If[total <= 0, Return[{1, 0}]];
    n = Length[row];
    pk = First@Ordering[row, -1];
    lo = pk; hi = pk; acc = row[[pk]];
    While[acc < windowAlpha total,
      lv = If[lo - 1 >= 1, row[[lo - 1]], -Infinity];
      rv = If[hi + 1 <= n, row[[hi + 1]], -Infinity];
      If[lv === -Infinity && rv === -Infinity, Break[]];
      If[lv >= rv, lo--; acc += row[[lo]], hi++; acc += row[[hi]]]];
    {lo, hi}
  ];
  If[windowAlpha > 0,
    matrix = (Module[{w = centralWindow[#]},
       Table[If[w[[1]] <= j <= w[[2]], #[[j]], 0.], {j, nIntervals}]] & /@ matrix)];

  (* ---- trim to the union of the per-row windows ----
     Keep columns where ANY row is nonzero; drop the all-zero margins. v_min
     rows are trimmed identically so the columns stay aligned. *)
  firstCol = SelectFirst[Range[nIntervals], Total[Abs[matrix[[All, #]]]] != 0 &, 1];
  lastCol  = SelectFirst[Range[nIntervals, 1, -1], Total[Abs[matrix[[All, #]]]] != 0 &, nIntervals];
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
