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
(*  Per-bin high-v_min cut, from a central window computed on the SMOOTH         *)
(*  response R(v_min) BEFORE integration: grow a window around the peak until it *)
(*  covers a fraction `alpha` of the response's area, then cut ONLY at its upper *)
(*  edge b -- keep [threshold, b] and remove only the decaying high-v_min tail   *)
(*  (the low rise is kept; R is zero below threshold anyway). The cut point is a *)
(*  property of the response, not of the v_min binning, so it is found on a fine *)
(*  grid. The matrix is then trimmed to the union of the per-bin windows.        *)
(*                                                                            *)
(*  This stage is light: it only imports the .wdx (self-contained) -- it does   *)
(*  NOT load the 01-06 pipeline.                                               *)
(*                                                                            *)
(*  Usage:                                                                     *)
(*  wolframscript -file 10_response_matrix.wl <name> <vminLo> <vminHi> <N> [alpha] [vLowFloor] *)
(*                                                                            *)
(*    <name>   base name of a .wdx in output/TES/response_functions/, or ALL   *)
(*    <vminLo> <vminHi>  integration range [km/s] (clipped to the function's   *)
(*                       domain if it exceeds it)                              *)
(*    <N>      number of equal v_min intervals (matrix columns)               *)
(*    [alpha]  central area fraction whose UPPER edge sets the high-v_min cut   *)
(*             (default 0.5). alpha = 0 disables the cut.                       *)
(*    [vLowFloor]  low-v_min cut [km/s]. Default (omitted) = the analytic        *)
(*             kinematic threshold sqrt(2 omega_thr / m_chi) of each config      *)
(*             (omega_thr = lowest bin edge, m_chi from the mass tag). A numeric  *)
(*             value overrides it for all configs.                              *)
(*                                                                            *)
(*  e.g.  wolframscript -file 10_response_matrix.wl ALL 1 800 1000 0.3          *)
(*        wolframscript -file 10_response_matrix.wl Al_q0M3_R5 1 800 1000 0.3 14 *)
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
(* Low-v_min floor [km/s]. Default (no 6th arg) = Automatic: each config is cut
   at its analytic kinematic threshold v_min = sqrt(2 omega_thr / m_chi), with
   omega_thr the lowest energy bin's lower edge and m_chi the DM mass (computed
   per config in buildMatrix). A numeric 6th arg overrides it for all configs. *)
vLowFloorCLI = If[Length[commandLineArgs] >= 6, ToExpression[commandLineArgs[[6]]], Automatic];
If[!(NumericQ[vminLoReq] && NumericQ[vminHiReq] && IntegerQ[nIntervals] && nIntervals >= 1
     && vminLoReq < vminHiReq && NumericQ[windowAlpha] && 0 <= windowAlpha < 1
     && (vLowFloorCLI === Automatic || (NumericQ[vLowFloorCLI] && vLowFloorCLI >= 0))),
  Print["error: need numeric vminLo < vminHi, integer N >= 1, 0 <= alpha < 1, vLowFloor >= 0"];
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
   vminLo, vminHi, edges, mids, matrix, vminRows, firstCol, stem, massTag, outDir,
   omegaThr, dmMassNat, vThr, floor},

  wdxPath = FileNameJoin[{functionDir, wdxBaseName <> ".wdx"}];
  If[!FileExistsQ[wdxPath], Print["skip (not found): ", wdxBaseName <> ".wdx"]; Return[]];

  data      = Import[wdxPath];
  responses = data["responses"];
  labels    = Keys[responses];
  fns       = Values[responses];
  nBins     = Length[fns];

  (* Low-v_min floor: analytic kinematic threshold v_min = sqrt(2 omega_thr/m_chi)
     unless the CLI overrides it. omega_thr = lowest bin's lower edge [eV], m_chi
     from the mass tag. The response is ~0 below this; the cut removes the
     resolution-smeared sub-threshold tail. *)
  omegaThr  = Min[First /@ (ToExpression /@ StringCases[#, NumberString] & /@ labels)] eV;
  dmMassNat = Switch[extractMass[wdxBaseName], "M1", 10 MeV, "M2", 100 MeV, "M3", 1000 MeV, _, $Failed];
  vThr      = If[NumericQ[dmMassNat] && omegaThr > 0, Sqrt[2 omegaThr/dmMassNat]/kps, 0.];
  floor     = If[vLowFloorCLI === Automatic, vThr, vLowFloorCLI];

  (* clip the requested range to the interpolation domain *)
  {domLo, domHi} = First[fns]["Domain"][[1]];
  vminLo = Max[vminLoReq, domLo];
  vminHi = Min[vminHiReq, domHi];
  If[{vminLo, vminHi} =!= {vminLoReq, vminHiReq},
    Print["  note: range clipped to function domain {", domLo, ", ", domHi, "}"]];

  edges = Subdivide[vminLo, vminHi, nIntervals];          (* N+1 edges *)
  mids  = (Most[edges] + Rest[edges]) / 2;                 (* N midpoints *)

  (* ---- per-bin central window on the SMOOTH response (before integration) ----
     The highest-density window {a,b} covering fraction windowAlpha of a
     response's area is a property of the smooth response function R(v_min), NOT
     of the v_min binning, so it is found here on a fine grid (independent of the
     matrix's N): grow outward from the peak, annexing the larger neighbour,
     until the kept area reaches windowAlpha of the total. windowAlpha = 0 keeps
     the full domain. *)
  nFine  = 4000;
  fineV  = Subdivide[vminLo, vminHi, nFine];
  fineDv = (vminHi - vminLo) / nFine;
  windowEdges[f_] := Module[{rs, total, m, pk, lo, hi, acc, lv, rv},
    rs = f /@ fineV;
    total = Total[rs] fineDv;
    If[total <= 0, Return[{vminHi, vminLo}]];            (* empty -> a > b *)
    m = Length[fineV];
    pk = First@Ordering[rs, -1];
    lo = pk; hi = pk; acc = rs[[pk]] fineDv;
    While[acc < windowAlpha total,
      lv = If[lo - 1 >= 1, rs[[lo - 1]], -Infinity];
      rv = If[hi + 1 <= m, rs[[hi + 1]], -Infinity];
      If[lv === -Infinity && rv === -Infinity, Break[]];
      If[lv >= rv, lo--; acc += rs[[lo]] fineDv, hi++; acc += rs[[hi]] fineDv]];
    {fineV[[lo]], fineV[[hi]]}
  ];
  winEdges = If[windowAlpha > 0, windowEdges /@ fns,
                ConstantArray[{vminLo, vminHi}, nBins]];

  (* M[i, j] = integral of bin-i response over interval j, with ONLY the high-
     v_min side cut at the window's upper edge b = winEdges[[i,2]]. The low side
     is kept down to the kinematic threshold (R is zero below it), i.e. keep
     [threshold, b] and remove only the decaying high-v_min tail. This avoids the
     low-edge spike a two-sided cut can introduce. kps converts dv to natural. *)
  matrix   = kps Table[
    Module[{loE = Max[edges[[j]], floor],
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
