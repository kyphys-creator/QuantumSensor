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


(* ============================ Unit conversion ============================ *)
(* The kernels are computed in natural units (all quantities in GeV, see
   01_setup). The response d(curly R)/d omega' carries dimension GeV^-2, i.e.
   (mass * energy)^-1. Multiplying by the natural-unit values of kg and eV
   (i.e. dividing by the physical unit kg^-1 eV^-1 = 1/(kg eV)) expresses it
   numerically in kg^-1 eV^-1. No hard-coded factor: kg and eV come from
   01_setup, so the conversion stays consistent with the rest of the code.   *)
toPhysicalKgEv = kg eV;


(* ============================ Shared plot styling ============================ *)

(* One colour per DM mass (Wolfram colour scheme 116). *)
colorMchi10MeV  = ColorData[116, 1];
colorMchi100MeV = ColorData[116, 2];
colorMchi1GeV   = ColorData[116, 3];

(* Common axis labels. *)
axisLabelVmin = Style["\!\(\*SubscriptBox[\(v\), \(min\)]\) [km/s]", 20, SingleLetterItalics -> False];
axisLabelRate = Style[" d\[ScriptCapitalR]/d\[Omega]' [\!\(\*SuperscriptBox[\(kg\), \(-1\)]\) \!\(\*SuperscriptBox[\(eV\), \(-1\)]\)]", 20];

(* kernelCurveStyle[massColour] -> options for the data curves ONLY (no frame).
   Solid style = first curve, dashed = second. The frame, axis labels and title
   are added separately (by `framed` for single plots, or frameDecor at the
   combined Show). Keeping the curve plots frame-less is what prevents a
   spurious vertical line: when several plots are overlaid with Show, each
   sub-plot would otherwise draw its own axes, and on the log v_min axis a
   sub-plot's y-axis is placed at the left edge of *its* range. The 10 MeV
   plot begins near v_min ~ 10, so its y-axis appeared as a stray vertical
   line at v_min ~ 10. Axes -> False removes it; the single frame from
   frameDecor provides the box and ticks instead.                             *)
kernelCurveStyle[massColor_] := {
  MaxRecursion -> 0, PlotPoints -> 800,
  ScalingFunctions -> {"Log"},
  PlotStyle -> {{Thick, massColor}, {Dashed, massColor}},
  Axes -> False,
  PlotRange -> All
};

(* frameDecor[titleText] -> frame / axis-label / title options, applied once at
   the Show level so the whole figure has a single frame. *)
frameDecor[titleText_] := {
  Frame -> True,
  LabelStyle -> Directive[FontFamily -> "Times", Black, 20, Bold],
  FrameStyle -> {Black},
  FrameLabel -> {axisLabelVmin, axisLabelRate, Style[titleText, 20]},
  ImageSize -> Large
};

(* framed[plot, title] -> one kernel curve-plot decorated with frame + labels. *)
framed[plot_, titleText_] := Show[plot, Sequence @@ frameDecor[titleText]];


(* ============================ Heavy mediator (n = 0) ============================ *)
(* KerRAll = left kinematic branch, KerRAlr = right branch, evaluated over
   v_min on a log y-axis. Recoil energy E_R = 1 eV.                           *)

plotAlHeavy10MeV = Plot[
  toPhysicalKgEv {-(KerRAll[10 MeV][0][1 eV, TESsig][vmin kps]) +
     (KerRAlr[10 MeV][0][1 eV, TESsig][vmin kps])},
  {vmin, 0, 10000},
  Evaluate[Sequence @@ kernelCurveStyle[colorMchi10MeV]]];

plotAlHeavy100MeV = Plot[
  toPhysicalKgEv {KerRAll[100 MeV][0][1, TESsig*1][vmin kps],
   KerRAlr[100 MeV][0][1, TESsig*1][vmin kps]},
  {vmin, 0, 1000},
  Evaluate[Sequence @@ kernelCurveStyle[colorMchi100MeV]]];

plotAlHeavy1GeV = Plot[
  toPhysicalKgEv {KerRAll[1000 MeV][0][1, TESsig*1][vmin kps],
   KerRAlr[1000 MeV][0][1, TESsig*1][vmin kps]},
  {vmin, 0, 1000},
  Evaluate[Sequence @@ kernelCurveStyle[colorMchi1GeV]]];


(* ============================ Light mediator (n = 2) ============================ *)
(* Two curves per plot: left+right branch summed, at E_R = 0.1 eV (solid) and
   E_R = 1 eV (dashed).                                                       *)

plotAlLight10MeV = Plot[
  toPhysicalKgEv {(KerRAll[10 MeV][2][0.1, TESsig*.1][vmin kps] +
      KerRAlr[10 MeV][2][0.1, TESsig*.1][vmin kps]),
   (KerRAll[10 MeV][2][1, TESsig*1][vmin kps] +
      KerRAlr[10 MeV][2][1, TESsig*1][vmin kps])},
  {vmin, 0, 10000},
  Evaluate[Sequence @@ kernelCurveStyle[colorMchi10MeV]]];

plotAlLight100MeV = Plot[
  toPhysicalKgEv {(KerRAll[100 MeV][2][0.1, TESsig*.1][vmin kps] +
      KerRAlr[100 MeV][2][0.1, TESsig*.1][vmin kps]),
   (KerRAll[100 MeV][2][1, TESsig*1][vmin kps] +
      KerRAlr[100 MeV][2][1, TESsig*1][vmin kps])},
  {vmin, 0, 1000},
  Evaluate[Sequence @@ kernelCurveStyle[colorMchi100MeV]]];

plotAlLight1GeV = Plot[
  toPhysicalKgEv {(KerRAll[1000 MeV][2][0.1, TESsig*.1][vmin kps] +
      KerRAlr[1000 MeV][2][0.1, TESsig*.1][vmin kps]),
   (KerRAll[1000 MeV][2][1, TESsig*1][vmin kps] +
      KerRAlr[1000 MeV][2][1, TESsig*1][vmin kps])},
  {vmin, 0, 800},
  Evaluate[Sequence @@ kernelCurveStyle[colorMchi1GeV]]];


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
  Sequence @@ frameDecor["TES, heavy mediator"],
  Epilog -> annotationLabels, PlotRange -> All];

figAlLight = Show[
  plotAlLight10MeV, plotAlLight100MeV, plotAlLight1GeV,
  Sequence @@ frameDecor["TES, light mediator"],
  Epilog -> annotationLabels, PlotRange -> All];


(* ============================ Export every figure ============================ *)

allFigures = {
  "Al_heavy_10MeV.pdf"   -> framed[plotAlHeavy10MeV,
    "TES, heavy mediator, m\[Chi]=10MeV, \!\(\*SubscriptBox[\(E\), \(R\)]\)=1eV"],
  "Al_heavy_100MeV.pdf"  -> framed[plotAlHeavy100MeV,
    "TES, heavy mediator, m\[Chi]=100MeV, \!\(\*SubscriptBox[\(E\), \(R\)]\)=1eV"],
  "Al_heavy_1GeV.pdf"    -> framed[plotAlHeavy1GeV,
    "TES, heavy mediator, m\[Chi]=1GeV, \!\(\*SubscriptBox[\(E\), \(R\)]\)=1eV"],
  "Al_light_10MeV.pdf"   -> framed[plotAlLight10MeV,  "TES, light mediator"],
  "Al_light_100MeV.pdf"  -> framed[plotAlLight100MeV, "TES, light mediator"],
  "Al_light_1GeV.pdf"    -> framed[plotAlLight1GeV,   "TES, light mediator"],
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
