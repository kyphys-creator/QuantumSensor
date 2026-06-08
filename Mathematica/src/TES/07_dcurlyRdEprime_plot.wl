(* ========================================================================== *)
(*  TES - d(curly R)/dE' kernel plots                                         *)
(*  Migrated from 07_dcurlyRdEprime_plot.nb.                                   *)
(*                                                                            *)
(*  Plots the per-energy response kernel  d curly R / d omega'  versus the    *)
(*  minimum DM velocity v_min, for Al (Mermin), three DM masses               *)
(*  (10 MeV / 100 MeV / 1 GeV) and two mediator types:                        *)
(*     heavy mediator -> FDM index n = 0                                       *)
(*     light mediator -> FDM index n = 2                                       *)
(*                                                                            *)
(*  Every figure (6 individual + 2 combined) is written as a PDF to           *)
(*     output/TES/dcurlyRdEprime/ .                                            *)
(*                                                                            *)
(*  Naming convention:                                                        *)
(*     plotAl<Heavy|Light><mass>   single kernel plot                          *)
(*     figAl<Heavy|Light>          combined figure (3 masses overlaid)         *)
(*                                                                            *)
(*  Notes on the migration (behaviour-preserving):                            *)
(*   - The .nb defined the heavy-mediator plots twice (an E_R = 0.1 eV draft   *)
(*     immediately overwritten by an E_R = 1 eV version). Only the surviving   *)
(*     (E_R = 1 eV) definitions are kept, matching the notebook's final state. *)
(*   - The .nb also carried TiN plot/export cells, but the TiN kernels are     *)
(*     never defined in this TES notebook. They are kept commented out at the  *)
(*     bottom so nothing errors at run time.                                   *)
(*   - Original Export used FileNameJoin[Direc, ...] (missing braces) and      *)
(*     wrote into input/TES; both are corrected here.                          *)
(* ========================================================================== *)

(* ============================ Setup ============================ *)

fileDir = DirectoryName[$InputFileName];
mathDir = ParentDirectory[ParentDirectory[fileDir]];
inputDir  = FileNameJoin[{mathDir, "input", "TES"}];
figureDir = FileNameJoin[{mathDir, "output", "TES", "dcurlyRdEprime"}];
If[!DirectoryQ[figureDir], CreateDirectory[figureDir, CreateIntermediateDirectories -> True]];

Get[FileNameJoin[{fileDir, "06_response_defs.wl"}]];
SetDirectory[inputDir];
Print["inputDir  = ", inputDir];
Print["figureDir = ", figureDir];


(* ============================ Shared plot styling ============================ *)

(* One colour per DM mass (Wolfram colour scheme 116). *)
colorMchi10MeV  = ColorData[116, 1];
colorMchi100MeV = ColorData[116, 2];
colorMchi1GeV   = ColorData[116, 3];

(* Common axis labels. *)
axisLabelVmin = Style["\!\(\*SubscriptBox[\(v\), \(min\)]\) [km/s]", 20, SingleLetterItalics -> False];
axisLabelRate = Style[" d\[ScriptCapitalR]/d\[Omega]' ( Scaled )", 20];

(* kernelPlotStyle[massColour, titleText] -> common Plot option list.
   The solid style is the first curve, the dashed style the second.          *)
kernelPlotStyle[massColor_, titleText_] := {
  MaxRecursion -> 0, PlotPoints -> 800,
  ScalingFunctions -> {"Log"},
  PlotStyle -> {{Thick, massColor}, {Dashed, massColor}},
  Frame -> True,
  LabelStyle -> Directive[FontFamily -> "Times", Black, 20, Bold],
  FrameStyle -> {Black},
  FrameLabel -> {axisLabelVmin, axisLabelRate, Style[titleText, 20]},
  PlotRange -> All, ImageSize -> Large
};


(* ============================ Heavy mediator (n = 0) ============================ *)
(* KerRAll = left kinematic branch, KerRAlr = right branch, evaluated over
   v_min on a log y-axis. Recoil energy E_R = 1 eV.                           *)

plotAlHeavy10MeV = Plot[
  {-(KerRAll[10 MeV][0][1 eV, TESsig][vmin kps]) +
     (KerRAlr[10 MeV][0][1 eV, TESsig][vmin kps])},
  {vmin, 0, 10000},
  Evaluate[Sequence @@ kernelPlotStyle[colorMchi10MeV,
    "TES, heavy mediator, m\[Chi]=10MeV, \!\(\*SubscriptBox[\(E\), \(R\)]\)=1eV"]]];

plotAlHeavy100MeV = Plot[
  {KerRAll[100 MeV][0][1, TESsig*1][vmin kps],
   KerRAlr[100 MeV][0][1, TESsig*1][vmin kps]},
  {vmin, 0, 1000},
  Evaluate[Sequence @@ kernelPlotStyle[colorMchi100MeV,
    "TES, heavy mediator, m\[Chi]=100MeV, \!\(\*SubscriptBox[\(E\), \(R\)]\)=1eV"]]];

plotAlHeavy1GeV = Plot[
  {KerRAll[1000 MeV][0][1, TESsig*1][vmin kps],
   KerRAlr[1000 MeV][0][1, TESsig*1][vmin kps]},
  {vmin, 0, 1000},
  Evaluate[Sequence @@ kernelPlotStyle[colorMchi1GeV,
    "TES, heavy mediator, m\[Chi]=1GeV, \!\(\*SubscriptBox[\(E\), \(R\)]\)=1eV"]]];


(* ============================ Light mediator (n = 2) ============================ *)
(* Two curves per plot: left+right branch summed, at E_R = 0.1 eV (solid) and
   E_R = 1 eV (dashed).                                                       *)

plotAlLight10MeV = Plot[
  {(KerRAll[10 MeV][2][0.1, TESsig*.1][vmin kps] +
      KerRAlr[10 MeV][2][0.1, TESsig*.1][vmin kps]),
   (KerRAll[10 MeV][2][1, TESsig*1][vmin kps] +
      KerRAlr[10 MeV][2][1, TESsig*1][vmin kps])},
  {vmin, 0, 10000},
  Evaluate[Sequence @@ kernelPlotStyle[colorMchi10MeV, "TES, light mediator"]]];

plotAlLight100MeV = Plot[
  {(KerRAll[100 MeV][2][0.1, TESsig*.1][vmin kps] +
      KerRAlr[100 MeV][2][0.1, TESsig*.1][vmin kps]),
   (KerRAll[100 MeV][2][1, TESsig*1][vmin kps] +
      KerRAlr[100 MeV][2][1, TESsig*1][vmin kps])},
  {vmin, 0, 1000},
  Evaluate[Sequence @@ kernelPlotStyle[colorMchi100MeV, "TES, light mediator"]]];

plotAlLight1GeV = Plot[
  {(KerRAll[1000 MeV][2][0.1, TESsig*.1][vmin kps] +
      KerRAlr[1000 MeV][2][0.1, TESsig*.1][vmin kps]),
   (KerRAll[1000 MeV][2][1, TESsig*1][vmin kps] +
      KerRAlr[1000 MeV][2][1, TESsig*1][vmin kps])},
  {vmin, 0, 800},
  Evaluate[Sequence @@ kernelPlotStyle[colorMchi1GeV, "TES, light mediator"]]];


(* ============================ Combined figures ============================ *)
(* In-figure text annotations (Inword from 02_functions_math) overlay the
   three mass curves. The labels use ImageScaled coordinates, so they float
   regardless of the data range.

   The original notebook also overlaid two horizontal reference lines at
   y = 0.175 / 0.095. Those only make sense for a normalised ("Scaled") plot:
   here the kernel values are raw (~10^-14), so lines at y ~ 0.1 sit ~13
   decades above the data and, on the shared log axis, stretched the y-limit
   so far that all three curves collapsed into an unreadable sliver. They are
   therefore dropped from the combined figures (see disabled block below).
   With only the kernel plots overlaid, the log y-axis auto-fits the data.    *)

annotationLabels = {
  Inword["\[Omega]'"][20][.16, .36][{Black, Bold}],
  Inword["0.1 eV"][18][.24, .3][Black],
  Inword["1 eV"][18][.228, .24][Black],
  Inword["\!\(\*SubscriptBox[\(m\), \(\[Chi]\)]\)"][20][.89, .85][Black],
  Inword["\!\(\*SuperscriptBox[\(10\), \(1\)]\) MeV"][18][.89, .78][colorMchi10MeV],
  Inword["\!\(\*SuperscriptBox[\(10\), \(2\)]\) MeV"][18][.89, .72][colorMchi100MeV],
  Inword["\!\(\*SuperscriptBox[\(10\), \(3\)]\) MeV"][18][.89, .66][colorMchi1GeV]
};

figAlHeavy = Show[
  plotAlHeavy10MeV, plotAlHeavy100MeV, plotAlHeavy1GeV,
  Epilog -> annotationLabels, PlotRange -> All];

figAlLight = Show[
  plotAlLight10MeV, plotAlLight100MeV, plotAlLight1GeV,
  Epilog -> annotationLabels, PlotRange -> All];


(* ============================ Export every figure ============================ *)

allFigures = {
  "Al_heavy_10MeV.pdf"   -> plotAlHeavy10MeV,
  "Al_heavy_100MeV.pdf"  -> plotAlHeavy100MeV,
  "Al_heavy_1GeV.pdf"    -> plotAlHeavy1GeV,
  "Al_light_10MeV.pdf"   -> plotAlLight10MeV,
  "Al_light_100MeV.pdf"  -> plotAlLight100MeV,
  "Al_light_1GeV.pdf"    -> plotAlLight1GeV,
  "Al_heavy_combined.pdf" -> figAlHeavy,
  "Al_light_combined.pdf" -> figAlLight
};

Do[
  Export[FileNameJoin[{figureDir, entry[[1]]}], entry[[2]]];
  Print["Saved: ", entry[[1]]],
  {entry, allFigures}];


(* ========================================================================== *)
(*  TiN section (disabled)                                                    *)
(*                                                                            *)
(*  The original notebook carried TiN equivalents, but the TiN kernel plots   *)
(*  are never defined in this TES notebook (they belong to the MKID pipeline). *)
(*  Kept for reference; enable only once the TiN kernel plots exist.           *)
(* --------------------------------------------------------------------------
   figTiNHeavy = Show[plotTiNHeavy10MeV, plotTiNHeavy100MeV, plotTiNHeavy1GeV,
     Epilog -> annotationLabels, PlotRange -> All];
   figTiNLight = Show[plotTiNLight10MeV, plotTiNLight100MeV, plotTiNLight1GeV,
     Epilog -> annotationLabels, PlotRange -> All];
   Export[FileNameJoin[{figureDir, "TiN_heavy_combined.pdf"}], figTiNHeavy];
   Export[FileNameJoin[{figureDir, "TiN_light_combined.pdf"}], figTiNLight];
   -------------------------------------------------------------------------- *)
