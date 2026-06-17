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
(*  Per-bin HIGH-v_min cut, computed on the SMOOTH response R(v_min) BEFORE      *)
(*  integration: keep [vminLo, b] where b is the OUTERMOST high-v_min at which R  *)
(*  >= a fraction `frac` of the peak value (default 0.2 -> cut the high-v tail    *)
(*  below 20% of the peak). Only the high side is cut; the low rise is kept (R is *)
(*  ~0 below the kinematic threshold, and all-zero low columns are trimmed). An   *)
(*  interior dip is kept, so a double-peaked high side is not split. b is a        *)
(*  property of the response, not of the v_min binning, so it is found on a fine   *)
(*  grid. The matrix is then trimmed to the union of the per-bin windows.          *)
(*                                                                            *)
(*  This stage is light: it only imports the .wdx (self-contained) -- it does   *)
(*  NOT load the 01-06 pipeline.                                               *)
(*                                                                            *)
(*  Usage:                                                                     *)
(*  wolframscript -file 10_response_matrix.wl <name> <vminLo> <vminHi> <N> [frac] *)
(*                                                                            *)
(*    <name>   base name of a .wdx in output/MKID/response_functions/, or ALL  *)
(*    <vminLo> <vminHi>  integration range [km/s] (clipped to the function's   *)
(*                       domain if it exceeds it)                              *)
(*    <N>      number of equal v_min intervals (matrix columns)               *)
(*    [frac]   peak-value fraction setting the HIGH-v_min cut: keep up to the   *)
(*             outermost v_min where R >= frac * peak (low side uncut)          *)
(*             (default 0.2). frac = 0 disables the cut.                        *)
(*                                                                            *)
(*  e.g.  wolframscript -file 10_response_matrix.wl ALL 1 800 1000 0.2          *)
(*        wolframscript -file 10_response_matrix.wl TiN_q0M3_R5 1 800 1000 0.2  *)
(*                                                                            *)
(*  Output (in output/MKID/response_matrix/<mass>/<name>_v<lo>-<hi>_N<N>/):    *)
(*    matrix.csv   pure numeric matrix (nBins x N)                             *)
(*    vmin.csv     per-column {v_low, v_high, v_mid}                           *)
(*    bins.csv     per-row bin energy labels                                    *)
(* ========================================================================== *)

(* ============================ Arguments ============================ *)

commandLineArgs = Rest[$ScriptCommandLine];
If[Length[commandLineArgs] < 4,
  Print["usage: wolframscript -file 10_response_matrix.wl <name|ALL> <vminLo> <vminHi> <N> [frac]"];
  Exit[1]];

nameArg    = commandLineArgs[[1]];
vminLoReq  = ToExpression[commandLineArgs[[2]]];
vminHiReq  = ToExpression[commandLineArgs[[3]]];
nIntervals = Round[ToExpression[commandLineArgs[[4]]]];
(* Per-bin high-v_min cut: keep [vminLo, b] where b is the outermost high-v_min
   at which the response stays >= fraction `peakFrac` of its peak value, cutting
   only the high-v_min tail below that level (default 0.2 -> cut where R has
   fallen to 20% of its peak; the low rise is kept). peakFrac = 0 disables it. *)
peakFrac = If[Length[commandLineArgs] >= 5, ToExpression[commandLineArgs[[5]]], 0.2];
If[!(NumericQ[vminLoReq] && NumericQ[vminHiReq] && IntegerQ[nIntervals] && nIntervals >= 1
     && vminLoReq < vminHiReq && NumericQ[peakFrac] && 0 <= peakFrac < 1),
  Print["error: need numeric vminLo < vminHi, integer N >= 1, 0 <= frac < 1"];
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

  (* ---- per-bin high-v_min cut on the SMOOTH response (before integration) ----
     Keep [vminLo, b] where b is the OUTERMOST high-v_min at which R >= peakFrac
     * peak. Only the high-v_min tail (below peakFrac of the peak) is cut; the low
     rise is kept (R ~ 0 below the kinematic threshold anyway, and all-zero low
     columns are trimmed below). b is a property of the smooth response, NOT of
     the v_min binning, so it is found on a fine grid. An interior dip is kept, so
     a double-peaked high side is not split. peakFrac = 0 keeps the full domain. *)
  nFine  = 4000;
  fineV  = Subdivide[vminLo, vminHi, nFine];
  windowEdges[f_] := Module[{rs, thr, idx},
    rs = f /@ fineV;
    If[Max[rs] <= 0, Return[{vminHi, vminLo}]];          (* empty -> a > b *)
    thr = peakFrac Max[rs];
    idx = Select[Range[Length[rs]], rs[[#]] >= thr &];
    {vminLo, fineV[[Last[idx]]]}      (* low: uncut; high: outermost crossing *)
  ];
  winEdges = If[peakFrac > 0, windowEdges /@ fns,
                ConstantArray[{vminLo, vminHi}, nBins]];

  (* M[i, j] = integral of bin-i response over interval j, clipped to the bin's
     window [vminLo, b] = winEdges[[i]]: only the high-v_min tail below
     peakFrac * peak is cut (b = winEdges[[i,2]]); the low rise is kept.
     kps converts dv to natural. *)
  matrix   = kps Table[
    Module[{loE = Max[edges[[j]], winEdges[[i, 1]]],
            hiE = Min[winEdges[[i, 2]], edges[[j + 1]]]},
      If[hiE > loE, integrateIF[fns[[i]], loE, hiE], 0.]],
    {i, nBins}, {j, nIntervals}];
  vminRows = N /@ Transpose[{Most[edges], Rest[edges], mids}];

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
