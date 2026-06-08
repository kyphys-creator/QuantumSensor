(* ::Package:: *)

(* ========================================================================== *)
(*  TES - Response-function plots                                              *)
(*  (formerly 09_response_matrix_plot; repurposed.)                           *)
(*                                                                            *)
(*  Plots the response functions that 08_response_functions.wl saved as .wdx.  *)
(*  For each .wdx (one DM mass + bin group) it overlays the per-energy-bin     *)
(*  InterpolatingFunctions  R_bin(v_min)  and writes a PDF.                    *)
(*                                                                            *)
(*  This stage is intentionally light: it does NOT load the 01-06 pipeline.    *)
(*  The .wdx files are self-contained (each holds the interpolating functions  *)
(*  plus metadata), so plotting is fast and decoupled from the heavy 08 run.   *)
(*                                                                            *)
(*  Usage:                                                                     *)
(*      wolframscript -file 09_response_function_plot.wl            (all files) *)
(*      wolframscript -file 09_response_function_plot.wl M3         (substring) *)
(*      wolframscript -file 09_response_function_plot.wl R5                     *)
(*  The optional argument keeps only .wdx whose name contains that substring.  *)
(*                                                                            *)
(*  Input:  output/TES/response_functions/*.wdx                                *)
(*  Output: output/TES/response_function_plots/*.pdf                           *)
(* ========================================================================== *)

(* ============================ Arguments ============================ *)

commandLineArgs = Rest[$ScriptCommandLine];
nameFilter = If[Length[commandLineArgs] >= 1, commandLineArgs[[1]], ""];


(* ============================ Setup ============================ *)

fileDir     = DirectoryName[$InputFileName];
mathDir     = ParentDirectory[ParentDirectory[fileDir]];
functionDir = FileNameJoin[{mathDir, "output", "TES", "response_functions"}];
plotDir     = FileNameJoin[{mathDir, "output", "TES", "response_function_plots"}];
If[!DirectoryQ[plotDir], CreateDirectory[plotDir, CreateIntermediateDirectories -> True]];
Print["functionDir = ", functionDir];
Print["plotDir     = ", plotDir];


(* ============================ Styling ============================ *)

axisLabelVmin = Style["\!\(\*SubscriptBox[\(v\), \(min\)]\) [km/s]", 16, SingleLetterItalics -> False];
axisLabelRate = Style["\!\(\*SubscriptBox[\(R\), \(bin\)]\)(\!\(\*SubscriptBox[\(v\), \(min\)]\))  [natural units]", 16];

(* One distinct colour per bin. *)
binColors[nBins_] := If[nBins <= 1, {ColorData["DarkRainbow"][0.5]},
  Table[ColorData["DarkRainbow"][(i - 1)/(nBins - 1)], {i, nBins}]];


(* ============================ Plot one file ============================ *)

(* plotResponseFile[wdxPath] -> overlay every bin's R_bin(v_min) and save a PDF. *)
plotResponseFile[wdxPath_] := Module[
  {data, responses, labels, fns, nBins, vlo, vhi, title, figure, outName},
  data      = Import[wdxPath];
  responses = data["responses"];
  labels    = Keys[responses];
  fns       = Values[responses];
  nBins     = Length[fns];
  {vlo, vhi} = data["vminRange"];
  title     = FileBaseName[wdxPath];

  figure = Plot[
    Evaluate[Through[fns[vmin]]],
    {vmin, vlo, vhi},
    PlotPoints -> 300, MaxRecursion -> 2,
    ScalingFunctions -> {"Log"},
    Axes -> False, Frame -> True,
    PlotStyle -> ({Thick, #} & /@ binColors[nBins]),
    PlotRange -> All,
    LabelStyle -> Directive[FontFamily -> "Times", Black, 16],
    FrameStyle -> Black,
    FrameLabel -> {axisLabelVmin, axisLabelRate, title},
    PlotLegends -> Placed[
      LineLegend[binColors[nBins], labels, LegendLabel -> "E' bin"], Right],
    ImageSize -> Large];

  outName = FileBaseName[wdxPath] <> ".pdf";
  Export[FileNameJoin[{plotDir, outName}], figure];
  Print["Saved: ", outName, "  (", nBins, " bins)"];
];


(* ============================ Dispatch ============================ *)

wdxFiles = FileNames["*.wdx", functionDir];
If[nameFilter =!= "",
  wdxFiles = Select[wdxFiles, StringContainsQ[FileBaseName[#], nameFilter] &]];

If[wdxFiles === {},
  Print["no .wdx files", If[nameFilter =!= "", " matching \"" <> nameFilter <> "\"", ""],
    " in ", functionDir],
  plotResponseFile /@ wdxFiles
];
