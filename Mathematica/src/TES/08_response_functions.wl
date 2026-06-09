(* ::Package:: *)

(* ========================================================================== *)
(*  TES - Binned velocity response (saved as interpolating functions)         *)
(*  (formerly 08_response_matrix; renamed since it now builds and saves        *)
(*  response functions, not a matrix.)                                         *)
(*                                                                            *)
(*  For each observed-energy bin [E1, E2] (in eV) this evaluates the binned    *)
(*  response  CRTES[mass, n][E1, E2][v_min]  on a v_min grid and builds an     *)
(*  InterpolatingFunction of it -- i.e. the response *as a function of v_min*.  *)
(*  The functions are saved to .wdx, so they can be re-imported later and      *)
(*  called at any v_min (no recomputation, no plot, no matrix).               *)
(*                                                                            *)
(*  Each .wdx holds one association:                                          *)
(*     <| "dmMass" -> ..., "fdmIndex" -> 0, "vminRange" -> {lo, hi},           *)
(*        "responses" -> <| "0.1-0.3 eV" -> InterpolatingFunction, ... |> |>   *)
(*                                                                            *)
(*  Reload & use, e.g.:                                                       *)
(*     data = Import["Al_q0M1_R5.wdx"];                                        *)
(*     f    = data["responses"]["0.1\[Dash]0.3 eV"];                          *)
(*     f[50.]      (* response at v_min = 50 km/s *)                           *)
(*                                                                            *)
(*  Run one group at a time (the curves are expensive), selected on the        *)
(*  command line:                                                             *)
(*                                                                            *)
(*      wolframscript -file 08_response_functions.wl bin5         (5 bins, all masses, heavy) *)
(*      wolframscript -file 08_response_functions.wl bin10        (10 bins, all masses, heavy)*)
(*      wolframscript -file 08_response_functions.wl bin5 M3      (5 bins, 1 GeV, heavy)       *)
(*      wolframscript -file 08_response_functions.wl bin5 ALL q2  (5 bins, all masses, light)  *)
(*                                                                            *)
(*  bin5  = 5 bins, 0.2 eV wide;  bin10 = 10 bins, 0.1 eV wide.                *)
(*  Optional 2nd arg M1 | M2 | M3 (default ALL) restricts to one mass.         *)
(*  Optional 3rd arg q0 | q2 (default q0) selects heavy (n=0) / light (n=2)    *)
(*  mediator. Al only. File name encodes the mediator: Al_q0M1_R5 / Al_q2M1_R5.*)
(*                                                                            *)
(*  Output: output/TES/response_functions/*.wdx                                *)
(* ========================================================================== *)

(* ============================ Arguments ============================ *)

commandLineArgs = Rest[$ScriptCommandLine];
binGroup     = If[Length[commandLineArgs] >= 1, ToLowerCase[commandLineArgs[[1]]], ""];
selectedMass = If[Length[commandLineArgs] >= 2, ToUpperCase[commandLineArgs[[2]]], "ALL"];
mediatorArg  = If[Length[commandLineArgs] >= 3, ToLowerCase[commandLineArgs[[3]]], "q0"];

(* Mediator -> {FDM index, file tag}. q0 = heavy (n=0), q2 = light (n=2). *)
{fdmIndex, qTag} = Switch[mediatorArg,
  "q0" | "heavy", {0, "q0"},
  "q2" | "light", {2, "q2"},
  _, {0, "q0"}];


(* ============================ Setup ============================ *)

fileDir     = DirectoryName[$InputFileName];
mathDir     = ParentDirectory[ParentDirectory[fileDir]];
inputDir    = FileNameJoin[{mathDir, "input", "TES"}];
functionDir = FileNameJoin[{mathDir, "output", "TES", "response_functions"}];
If[!DirectoryQ[functionDir], CreateDirectory[functionDir, CreateIntermediateDirectories -> True]];

Get[FileNameJoin[{fileDir, "06_response_defs.wl"}]];
SetDirectory[inputDir];
Print["functionDir = ", functionDir];

(* Note: parallelising the v_min samples was tried (ParallelMap with the defs
   distributed to subkernels) and was ~80x SLOWER -- CRTES closes over the large
   Mermin interpolation, so the per-call MathLink overhead dwarfs the ~8 ms of
   compute. It is therefore kept serial. The big win is the inlined gaussPDF in
   06 (see there), which already made CRTES ~3.6x faster. *)


(* ============================ Bins & masses ============================ *)

binEdges5  = Range[0.1, 1.1, 0.2];   (* {0.1,0.3,0.5,0.7,0.9,1.1} -> 5 bins  *)
binEdges10 = Range[0.1, 1.1, 0.1];   (* {0.1,0.2,...,1.1}         -> 10 bins *)

(* {DM mass, file tag, title text} for the heavy mediator. *)
dmMasses = {
  {10 MeV,   "M1", "m\[Chi]=10MeV"},
  {100 MeV,  "M2", "m\[Chi]=100MeV"},
  {1000 MeV, "M3", "m\[Chi]=1GeV"}
};


(* ============================ v_min sampling grid ============================ *)
(* The InterpolatingFunction is built on this grid [km/s]. The grid depends on
   the mediator: the light mediator (q2) response extends to higher v_min, so
   it uses a wider domain and more samples (same ~10 pts/km/s density).
   Increase nVminSamples for finer resolution (cost scales linearly).          *)

{vminMin, vminMax, nVminSamples} = If[fdmIndex == 2,
  {1, 2000, 20000},   (* light (q2) *)
  {1,  800,  8000}];  (* heavy (q0) *)
vminGrid = Subdivide[vminMin, vminMax, nVminSamples - 1];   (* nVminSamples points *)


(* ============================ Builders ============================ *)

(* Energy-range label per bin, e.g. "0.1-0.3 eV" -- used as the association key. *)
binEnergyLabels[binEdges_] := Table[
  ToString[binEdges[[i]]] <> "\[Dash]" <> ToString[binEdges[[i + 1]]] <> " eV",
  {i, Length[binEdges] - 1}];

(* buildResponseFunction[dmMass, fdmIndex][lowEdge, highEdge] -> an
   InterpolatingFunction of CRTES vs v_min for that energy bin. Order-1
   (piecewise-linear) interpolation is used deliberately: CRTES is >= 0
   everywhere, and a linear interpolant never leaves the [min, max] of its
   sample points, so it cannot overshoot or dip negative near the sharp kernel
   peaks the way order-2/3 would. With the dense v_min grid this is still
   smooth. *)
buildResponseFunction[dmMass_, fdmIndex_][lowEdge_, highEdge_] := Module[{values},
  values = CRTES[dmMass, fdmIndex][lowEdge, highEdge][#] & /@ vminGrid;
  Interpolation[Transpose[{vminGrid, values}], InterpolationOrder -> 1]
];


(* ============================ Run a group ============================ *)

runGroup[binTag_, binEdges_] := Do[
  Module[{dmMass, massTag, massTitle, nBins, responses, payload, fileName},
    {dmMass, massTag, massTitle} = massEntry;
    If[selectedMass === "ALL" || selectedMass === massTag,
      nBins = Length[binEdges] - 1;
      responses = AssociationThread[
        binEnergyLabels[binEdges] ->
        Table[
          buildResponseFunction[dmMass, fdmIndex][binEdges[[i]], binEdges[[i + 1]]],
          {i, nBins}]];
      payload = <|
        "dmMass"    -> dmMass,
        "fdmIndex"  -> fdmIndex,
        "vminRange" -> {vminMin, vminMax},
        "responses" -> responses |>;
      fileName = "Al_" <> qTag <> massTag <> "_" <> binTag <> ".wdx";
      Export[FileNameJoin[{functionDir, fileName}], payload];
      Print["Saved: ", fileName, "  (", nBins, " bins)"];
    ];
  ],
  {massEntry, dmMasses}
];


(* ============================ Dispatch ============================ *)

Which[
  binGroup === "bin5",  runGroup["R5",  binEdges5],
  binGroup === "bin10", runGroup["R10", binEdges10],
  True, Print[
    "usage: wolframscript -file 08_response_functions.wl [bin5|bin10] [M1|M2|M3|ALL] [q0|q2]"]
];
