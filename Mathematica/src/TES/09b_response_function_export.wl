(* ::Package:: *)

(* ========================================================================== *)
(*  TES - Response-function CSV export                                         *)
(*                                                                            *)
(*  Companion to 09_response_function_plot.wl: instead of a PDF, it writes a   *)
(*  plain CSV per .wdx so the Python side can plot the response functions over *)
(*  their FULL domain (the matrices/eta CSVs only carry the windowed grid).    *)
(*                                                                            *)
(*  For each .wdx (one DM mass + bin group) it samples every per-energy-bin    *)
(*  InterpolatingFunction R_bin(v_min) on a log-spaced v_min grid spanning the *)
(*  function's own domain, and writes                                          *)
(*                                                                            *)
(*      vmin, <bin1 label>, <bin2 label>, ...                                  *)
(*      v_1 , R_1*kg       , R_2*kg      , ...                                  *)
(*                                                                            *)
(*  The stored values are R_bin*kg, i.e. the [kg^-1] quantity that             *)
(*  09_response_function_plot.wl plots, so Python can plot the column as-is.    *)
(*                                                                            *)
(*  Usage:                                                                     *)
(*      wolframscript -file 09b_response_function_export.wl            (all)    *)
(*      wolframscript -file 09b_response_function_export.wl M3         (filter) *)
(*                                                                            *)
(*  Input:  output/TES/response_functions/*.wdx                                *)
(*  Output: output/TES/response_functions_csv/*.csv                            *)
(* ========================================================================== *)

commandLineArgs = Rest[$ScriptCommandLine];
nameFilter = If[Length[commandLineArgs] >= 1, commandLineArgs[[1]], ""];
nSamples   = If[Length[commandLineArgs] >= 2, Round[ToExpression[commandLineArgs[[2]]]], 600];

fileDir     = DirectoryName[$InputFileName];
mathDir     = ParentDirectory[ParentDirectory[fileDir]];
functionDir = FileNameJoin[{mathDir, "output", "TES", "response_functions"}];
outDir      = FileNameJoin[{mathDir, "output", "TES", "response_functions_csv"}];
If[!DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];

Get[FileNameJoin[{fileDir, "01_setup.wl"}]];   (* defines kg *)

Print["functionDir = ", functionDir];
Print["outDir      = ", outDir];

exportFile[wdxPath_] := Module[
  {data, responses, labels, fns, dom, lo, hi, vs, table, header, outName},
  data      = Import[wdxPath];
  responses = data["responses"];
  labels    = Keys[responses];
  fns       = Values[responses];
  {lo, hi}  = First[fns]["Domain"][[1]];
  (* log-spaced sample grid over the full interpolation domain *)
  vs    = Exp[Subdivide[Log[N@lo], Log[N@hi], nSamples - 1]];
  (* each row: {v, R_bin1(v)*kg, R_bin2(v)*kg, ...}  (the [kg^-1] plot quantity) *)
  table = Table[Prepend[(#[v] & /@ fns) kg, v], {v, vs}];
  header  = Prepend[labels, "vmin"];
  outName = FileBaseName[wdxPath] <> ".csv";
  Export[FileNameJoin[{outDir, outName}], Prepend[table, header]];
  Print["Saved: ", outName, "  (", Length[fns], " bins x ", nSamples, " pts, domain ",
        N@lo, "-", N@hi, ")"];
];

wdxFiles = FileNames["*.wdx", functionDir];
If[nameFilter =!= "",
  wdxFiles = Select[wdxFiles, StringContainsQ[FileBaseName[#], nameFilter] &]];

If[wdxFiles === {},
  Print["no .wdx files", If[nameFilter =!= "", " matching \"" <> nameFilter <> "\"", ""],
    " in ", functionDir],
  exportFile /@ wdxFiles
];
