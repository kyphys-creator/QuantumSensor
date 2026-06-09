(* ::Package:: *)

(* ========================================================================== *)
(*  MKID - d(curly R)/dE' kernel plots                                        *)
(*  TiN counterpart of TES/07_dcurlyRdEprime_plot.wl.                          *)
(*                                                                            *)
(*  Plots the per-energy response kernel  d curly R / d omega'  versus the    *)
(*  minimum DM velocity v_min, for TiN (analytic Lindhard), three DM masses   *)
(*  (10 MeV / 100 MeV / 1 GeV) and two mediator types:                        *)
(*     heavy mediator -> FDM index n = 0                                       *)
(*     light mediator -> FDM index n = 2                                       *)
(*                                                                            *)
(*  Every figure (6 individual + 2 combined) is written as a PDF to           *)
(*     output/MKID/dcurlyRdEprime/ .                                           *)
(*                                                                            *)
(*  Naming convention:                                                        *)
(*     plotTiN<Heavy|Light><mass>   single kernel plot                         *)
(*     figTiN<Heavy|Light>          combined figure (3 masses overlaid)        *)
(*                                                                            *)
(*  The TiN kernels (KerRTiNl / KerRTiNr) come from 06_response_defs.wl;       *)
(*  the resolution MKIDsig and units (kg, eV, kps, ...) come from 01/05.       *)
(* ========================================================================== *)

(* ============================ Setup ============================ *)

fileDir = DirectoryName[$InputFileName];
mathDir = ParentDirectory[ParentDirectory[fileDir]];
inputDir  = FileNameJoin[{mathDir, "input", "MKID"}];
figureDir = FileNameJoin[{mathDir, "output", "MKID", "dcurlyRdEprime"}];
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
   See TES/07 for the rationale behind Axes -> False (it suppresses a spurious
   vertical line when several log-axis sub-plots are overlaid with Show; the
   single frame from frameDecor provides the box and ticks instead).          *)
kernelCurveStyle[massColor_] := {
  MaxRecursion -> 0, PlotPoints -> 800,
  ScalingFunctions -> {"Log"},
  (* The two recoil energies share the mass colour but differ in opacity:
     E_R = 0.2 eV is the bold (opaque) curve, E_R = 1 eV the faded one. *)
  PlotStyle -> {{Thick, massColor}, {Thick, Opacity[0.4], massColor}},
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
(* KerRTiNl = left kinematic branch, KerRTiNr = right branch, evaluated over
   v_min on a log y-axis. Recoil energies E_R = 0.2 eV (bold) / 1 eV (faded).  *)

plotTiNHeavy10MeV = Plot[
  toPhysicalKgEv {(KerRTiNl[10 MeV][0][0.2 eV, MKIDsig][vmin kps]) +
     (KerRTiNr[10 MeV][0][0.2 eV, MKIDsig][vmin kps]),
     (KerRTiNl[10 MeV][0][1 eV, MKIDsig][vmin kps]) +
     (KerRTiNr[10 MeV][0][1 eV, MKIDsig][vmin kps])},
  {vmin, 1, 10000},
  Evaluate[Sequence @@ kernelCurveStyle[colorMchi10MeV]]];

plotTiNHeavy100MeV = Plot[
  toPhysicalKgEv {(KerRTiNl[100 MeV][0][0.2 eV, MKIDsig][vmin kps]) +
     (KerRTiNr[100 MeV][0][0.2 eV, MKIDsig][vmin kps]),
     (KerRTiNl[100 MeV][0][1 eV, MKIDsig][vmin kps]) +
     (KerRTiNr[100 MeV][0][1 eV, MKIDsig][vmin kps])},
  {vmin, 1, 10000},
  Evaluate[Sequence @@ kernelCurveStyle[colorMchi100MeV]]];

plotTiNHeavy1GeV = Plot[
  toPhysicalKgEv {(KerRTiNl[1000 MeV][0][0.2 eV, MKIDsig][vmin kps]) +
     (KerRTiNr[1000 MeV][0][0.2 eV, MKIDsig][vmin kps]),
     (KerRTiNl[1000 MeV][0][1 eV, MKIDsig][vmin kps]) +
     (KerRTiNr[1000 MeV][0][1 eV, MKIDsig][vmin kps])},
  {vmin, 1, 10000},
  Evaluate[Sequence @@ kernelCurveStyle[colorMchi1GeV]]];


(* ============================ Light mediator (n = 2) ============================ *)
(* Two curves per plot: left+right branch summed, at E_R = 0.2 eV (bold) and
   E_R = 1 eV (faded).                                                        *)

plotTiNLight10MeV = Plot[
  toPhysicalKgEv {(KerRTiNl[10 MeV][2][0.2 eV, MKIDsig][vmin kps]) +
      (KerRTiNr[10 MeV][2][0.2 eV, MKIDsig][vmin kps]),
   (KerRTiNl[10 MeV][2][1 eV, MKIDsig][vmin kps]) +
      (KerRTiNr[10 MeV][2][1 eV, MKIDsig][vmin kps])},
  {vmin, 1, 10000},
  Evaluate[Sequence @@ kernelCurveStyle[colorMchi10MeV]]];

plotTiNLight100MeV = Plot[
  toPhysicalKgEv {(KerRTiNl[100 MeV][2][0.2 eV, MKIDsig][vmin kps]) +
      (KerRTiNr[100 MeV][2][0.2 eV, MKIDsig][vmin kps]),
   (KerRTiNl[100 MeV][2][1 eV, MKIDsig][vmin kps]) +
      (KerRTiNr[100 MeV][2][1 eV, MKIDsig][vmin kps])},
  {vmin, 1, 10000},
  Evaluate[Sequence @@ kernelCurveStyle[colorMchi100MeV]]];

plotTiNLight1GeV = Plot[
  toPhysicalKgEv {(KerRTiNl[1000 MeV][2][0.2 eV, MKIDsig][vmin kps]) +
      (KerRTiNr[1000 MeV][2][0.2 eV, MKIDsig][vmin kps]),
   (KerRTiNl[1000 MeV][2][1 eV, MKIDsig][vmin kps]) +
      (KerRTiNr[1000 MeV][2][1 eV, MKIDsig][vmin kps])},
  {vmin, 1, 10000},
  Evaluate[Sequence @@ kernelCurveStyle[colorMchi1GeV]]];


(* ============================ Combined figures ============================ *)
(* In-figure text annotations (Inword from 02_functions_math) overlay the
   three mass curves. The labels use ImageScaled coordinates, so they float
   regardless of the data range. (See TES/07 for why the legacy normalised
   reference lines are intentionally dropped here.)                           *)

(* legendEntry[lineStyle, label][x, y] -> a legend row combining a short line
   swatch and its text label, left-anchored at the ImageScaled position (x, y). *)
legendEntry[lineStyle_, label_][x_, y_] := Inset[
  Row[{
    Graphics[{Black, Thick, lineStyle, Line[{{0, 0}, {1, 0}}]},
      ImageSize -> 30, AspectRatio -> 0.16],
    Style[" " <> label, 18, Black, FontFamily -> "Times"]
  }],
  ImageScaled[{x, y}], {Left, Center}];

annotationLabels = {
  Inword["\[Omega]'"][20][.16, .36][{Black, Bold}],
  (* E_R legend: bold line + "0.2 eV", faded line + "1 eV" *)
  legendEntry[Opacity[1], "0.2 eV"][.18, .3],
  legendEntry[Opacity[0.4], "1 eV"][.18, .24],
  Inword["\!\(\*SubscriptBox[\(m\), \(\[Chi]\)]\)"][20][.89, .85][Black],
  Inword["\!\(\*SuperscriptBox[\(10\), \(1\)]\) MeV"][18][.89, .78][colorMchi10MeV],
  Inword["\!\(\*SuperscriptBox[\(10\), \(2\)]\) MeV"][18][.89, .72][colorMchi100MeV],
  Inword["\!\(\*SuperscriptBox[\(10\), \(3\)]\) MeV"][18][.89, .66][colorMchi1GeV]
};

figTiNHeavy = Show[
  plotTiNHeavy10MeV, plotTiNHeavy100MeV, plotTiNHeavy1GeV,
  Sequence @@ frameDecor["MKID, heavy mediator"],
  Epilog -> annotationLabels, PlotRange -> All];

figTiNLight = Show[
  plotTiNLight10MeV, plotTiNLight100MeV, plotTiNLight1GeV,
  Sequence @@ frameDecor["MKID, light mediator"],
  Epilog -> annotationLabels, PlotRange -> All];


(* ============================ Export every figure ============================ *)

allFigures = {
  "TiN_heavy_10MeV.pdf"   -> framed[plotTiNHeavy10MeV,
    "MKID, heavy mediator, m\[Chi]=10MeV, \!\(\*E'\)=0.2eV, 1eV"],
  "TiN_heavy_100MeV.pdf"  -> framed[plotTiNHeavy100MeV,
    "MKID, heavy mediator, m\[Chi]=100MeV, \!\(\*E'\)=0.2eV, 1eV"],
  "TiN_heavy_1GeV.pdf"    -> framed[plotTiNHeavy1GeV,
    "MKID, heavy mediator, m\[Chi]=1GeV, \!\(\*E'\)=0.2eV, 1eV"],
  "TiN_light_10MeV.pdf"   -> framed[plotTiNLight10MeV,  "MKID, light mediator"],
  "TiN_light_100MeV.pdf"  -> framed[plotTiNLight100MeV, "MKID, light mediator"],
  "TiN_light_1GeV.pdf"    -> framed[plotTiNLight1GeV,   "MKID, light mediator"],
  "TiN_heavy_combined.pdf" -> figTiNHeavy,
  "TiN_light_combined.pdf" -> figTiNLight
};

Do[
  Export[FileNameJoin[{figureDir, entry[[1]]}], entry[[2]]];
  Print["Saved: ", entry[[1]]],
  {entry, allFigures}];
