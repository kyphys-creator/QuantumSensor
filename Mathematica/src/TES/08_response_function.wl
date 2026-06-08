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
(*      wolframscript -file 08_response_function.wl bin5      (5 bins, all masses)   *)
(*      wolframscript -file 08_response_function.wl bin10     (10 bins, all masses)  *)
(*      wolframscript -file 08_response_function.wl bin5 M3   (5 bins, 1 GeV only)   *)
(*                                                                            *)
(*  bin5  = 5 bins, 0.2 eV wide;  bin10 = 10 bins, 0.1 eV wide.                *)
(*  Optional 2nd arg M1 | M2 | M3 restricts to one mass (10 MeV / 100 / 1 GeV).*)
(*  Heavy mediator (q0, fdmIndex=0) and Al only, as in the original; light      *)
(*  (fdmIndex=2) can be added by passing fdmIndex / extending dmMasses.         *)
(*                                                                            *)
(*  Output: output/TES/response_functions/*.wdx                                *)
(* ========================================================================== *)

(* ============================ Arguments ============================ *)

commandLineArgs = Rest[$ScriptCommandLine];
binGroup     = If[Length[commandLineArgs] >= 1, ToLowerCase[commandLineArgs[[1]]], ""];
selectedMass = If[Length[commandLineArgs] >= 2, ToUpperCase[commandLineArgs[[2]]], "ALL"];


(* ============================ Setup ============================ *)

fileDir     = DirectoryName[$InputFileName];
mathDir     = ParentDirectory[ParentDirectory[fileDir]];
inputDir    = FileNameJoin[{mathDir, "input", "TES"}];
functionDir = FileNameJoin[{mathDir, "output", "TES", "response_functions"}];
If[!DirectoryQ[functionDir], CreateDirectory[functionDir, CreateIntermediateDirectories -> True]];

Get[FileNameJoin[{fileDir, "06_response_defs.wl"}]];
SetDirectory[inputDir];
Print["functionDir = ", functionDir];


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
(* The InterpolatingFunction is built on this grid [km/s]. Increase nVminSamples
   for finer resolution (and longer run time); the cost scales linearly.       *)

vminMin = 1;  vminMax = 800;  nVminSamples = 800;
vminGrid = Subdivide[vminMin, vminMax, nVminSamples - 1];   (* nVminSamples points *)


(* ============================ Builders ============================ *)

(* Energy-range label per bin, e.g. "0.1-0.3 eV" -- used as the association key. *)
binEnergyLabels[binEdges_] := Table[
  ToString[binEdges[[i]]] <> "\[Dash]" <> ToString[binEdges[[i + 1]]] <> " eV",
  {i, Length[binEdges] - 1}];

(* buildResponseFunction[dmMass, fdmIndex][lowEdge, highEdge] -> an
   InterpolatingFunction of CRTES vs v_min for that energy bin. Order-1
   (piecewise-linear) interpolation, matching the rest of the codebase and
   robust to the sharp kernel features. *)
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
          buildResponseFunction[dmMass, 0][binEdges[[i]], binEdges[[i + 1]]],
          {i, nBins}]];
      payload = <|
        "dmMass"    -> dmMass,
        "fdmIndex"  -> 0,
        "vminRange" -> {vminMin, vminMax},
        "responses" -> responses |>;
      fileName = "Al_q0" <> massTag <> "_" <> binTag <> ".wdx";
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
    "usage: wolframscript -file 08_response_function.wl [bin5|bin10] [M1|M2|M3]"]
];
