(* ::Package:: *)

(* ========================================================================== *)
(*  MKID - Binned velocity response (saved as interpolating functions)        *)
(*  TiN counterpart of TES/08_response_functions.wl.                           *)
(*                                                                            *)
(*  For each observed-energy bin [E1, E2] (in eV) this evaluates the binned    *)
(*  response  CRTiN[mass, n][E1, E2][v_min]  on a v_min grid and builds an     *)
(*  InterpolatingFunction of it -- i.e. the response *as a function of v_min*.  *)
(*  The functions are saved to .wdx, so they can be re-imported later and      *)
(*  called at any v_min (no recomputation, no plot, no matrix).               *)
(*                                                                            *)
(*  Each .wdx holds one association:                                          *)
(*     <| "dmMass" -> ..., "fdmIndex" -> 0, "vminRange" -> {lo, hi},           *)
(*        "responses" -> <| "0.1-0.3 eV" -> InterpolatingFunction, ... |> |>   *)
(*                                                                            *)
(*  Reload & use, e.g.:                                                       *)
(*     data = Import["TiN_q0M1_R5.wdx"];                                       *)
(*     f    = data["responses"]["0.1\[Dash]0.3 eV"];                          *)
(*     f[50.]      (* response at v_min = 50 km/s *)                           *)
(*                                                                            *)
(*  Run one group at a time (the curves are expensive), selected on the        *)
(*  command line:                                                             *)
(*                                                                            *)
(*      wolframscript -file 08_response_functions.wl bin5         (5 bins, all masses, heavy) *)
(*      wolframscript -file 08_response_functions.wl bin9         (9 bins, all masses, heavy) *)
(*      wolframscript -file 08_response_functions.wl bin5 M3      (5 bins, 1 GeV, heavy)       *)
(*      wolframscript -file 08_response_functions.wl bin5 ALL q2  (5 bins, all masses, light)  *)
(*                                                                            *)
(*  MKID threshold = 0.2 eV, so bins start at 0.2 eV:                          *)
(*    bin5  -> [0.2,0.3],[0.3,0.5],[0.5,0.7],[0.7,0.9],[0.9,1.1]  (5 bins)      *)
(*    bin9  -> [0.2,0.3],[0.3,0.4],...,[1.0,1.1]                  (9 bins)      *)
(*  Optional 2nd arg M1 | M2 | M3 (default ALL) restricts to one mass.         *)
(*  Optional 3rd arg q0 | q2 (default q0) selects heavy (n=0) / light (n=2)    *)
(*  mediator. TiN only. File name encodes the mediator: TiN_q0M1_R5 / TiN_q2M1_R5.*)
(*                                                                            *)
(*  Output: output/MKID/response_functions/*.wdx                               *)
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
inputDir    = FileNameJoin[{mathDir, "input", "MKID"}];
functionDir = FileNameJoin[{mathDir, "output", "MKID", "response_functions"}];
If[!DirectoryQ[functionDir], CreateDirectory[functionDir, CreateIntermediateDirectories -> True]];

Get[FileNameJoin[{fileDir, "06_response_defs.wl"}]];
SetDirectory[inputDir];
Print["functionDir = ", functionDir];

(* Kept serial deliberately (see TES/08): the per-call MathLink overhead of
   distributing CRTiN to subkernels dwarfs the compute, so ParallelMap is a
   large net loss. The inlined gaussPDF in 06 already speeds CRTiN up ~4x.    *)


(* ============================ Bins & masses ============================ *)

(* MKID detection threshold is 0.2 eV, so every bin starts at or above 0.2 eV
   (unlike TES, which goes down to 0.1 eV). The lowest bin is a narrow 0.1 eV
   bin [0.2, 0.3]; the rest of bin5 is 0.2 eV wide. *)
binEdges5  = {0.2, 0.3, 0.5, 0.7, 0.9, 1.1};   (* [0.2,0.3],[0.3,0.5],[0.5,0.7],[0.7,0.9],[0.9,1.1] -> 5 bins *)
binEdges9  = Range[0.2, 1.1, 0.1];             (* {0.2,0.3,...,1.1} -> 9 bins of 0.1 eV (tag R9)         *)

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
   InterpolatingFunction of CRTiN vs v_min for that energy bin. Order-1
   (piecewise-linear) interpolation is used deliberately: CRTiN is >= 0
   everywhere, and a linear interpolant never leaves the [min, max] of its
   sample points, so it cannot overshoot or dip negative near the sharp kernel
   peaks the way order-2/3 would. With the dense v_min grid this is still
   smooth. *)
buildResponseFunction[dmMass_, fdmIndex_][lowEdge_, highEdge_] := Module[{values},
  values = CRTiN[dmMass, fdmIndex][lowEdge, highEdge][#] & /@ vminGrid;
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
      fileName = "TiN_" <> qTag <> massTag <> "_" <> binTag <> ".wdx";
      Export[FileNameJoin[{functionDir, fileName}], payload];
      Print["Saved: ", fileName, "  (", nBins, " bins)"];
    ];
  ],
  {massEntry, dmMasses}
];


(* ============================ Dispatch ============================ *)

Which[
  binGroup === "bin5", runGroup["R5", binEdges5],
  binGroup === "bin9", runGroup["R9", binEdges9],
  True, Print[
    "usage: wolframscript -file 08_response_functions.wl [bin5|bin9] [M1|M2|M3|ALL] [q0|q2]"]
];
