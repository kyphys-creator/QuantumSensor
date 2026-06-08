#!/usr/bin/env wolframscript

(* ==========================================================================
   test_definitions.wls — Smoke test for the 01-06 .wl pipeline.

   Loads <det>/06_response_defs.wl (which transitively Gets 01-05) and checks:
     (A) the load itself raises no messages / errors,
     (B) every expected symbol from each stage is actually defined,
     (C) representative calls return finite real numbers (numerical smoke).

   Usage:
     wolframscript -file test_definitions.wls [TES|MKID]

   Default detector is TES. A human-readable report and a machine-readable
   association are written to  output/<det>/test/.
   Exit code is 0 only if every check passes.
   ========================================================================== *)

(* ---------------------------------------------------------------- detector *)
args = Rest[$ScriptCommandLine];
det  = If[Length[args] >= 1, ToUpperCase[args[[1]]], "TES"];
If[!MemberQ[{"TES", "MKID"}, det],
  Print["unknown detector: ", det, " (expected TES or MKID)"]; Exit[2]];

(* ------------------------------------------------------------------- paths *)
scriptDir = DirectoryName[$InputFileName];                 (* .../src/tests   *)
srcDir    = ParentDirectory[scriptDir];                    (* .../src         *)
mathDir   = ParentDirectory[srcDir];                       (* .../Mathematica *)
detSrc    = FileNameJoin[{srcDir, det, "06_response_defs.wl"}];
outDir    = FileNameJoin[{mathDir, "output", det, "test"}];
If[!DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];

Print["==== test_definitions :: ", det, " ===="];
Print["loading: ", detSrc];

(* ---------------------------------------------------------------- recorder *)
results = {};                       (* list of <|name, status, detail|>      *)
record[name_, status_, detail_:""] := (
  AppendTo[results, <|"name" -> name, "status" -> status, "detail" -> detail|>];
  Print["  [", status, "] ", name, If[detail === "", "", "  -- " <> ToString[detail]]];
);

(* ------------------------------------------------------------- (A) loading *)
(* Capture any messages emitted while loading, but keep going. *)
loadOK =
  Check[
    Quiet[Get[detSrc]; True] // TrueQ,
    False
  ];
record["load/06_response_defs",
  If[loadOK, "PASS", "FAIL"],
  If[loadOK, detSrc, "Get raised an error — see kernel output above"]];

(* If the load outright failed we cannot check anything else meaningfully. *)
If[!loadOK,
  Export[FileNameJoin[{outDir, "report.txt"}],
    "LOAD FAILED for " <> det <> " — " <> detSrc, "Text"];
  Print["LOAD FAILED — aborting."]; Exit[1]];

(* ------------------------------------------------ defined? symbol-level check *)
(* A symbol counts as "defined" if it carries Own/Down/Sub values or is a
   bound numeric constant. ToExpression keeps the name held in Global`.        *)
SetAttributes[definedQ, HoldFirst];
definedQ[name_String] := ToExpression[name, InputForm,
  Function[s,
    (OwnValues[s]  =!= {}) ||
    (DownValues[s] =!= {}) ||
    (SubValues[s]  =!= {}) ||
    NumericQ[Quiet[s]],
    HoldAll]];

checkDefined[stage_String, names_List] :=
  Module[{missing},
    missing = Select[names, !definedQ[#] &];
    record["defs/" <> stage,
      If[missing === {}, "PASS", "FAIL"],
      If[missing === {}, ToString[Length[names]] <> " symbols",
        "missing: " <> StringRiffle[missing, ", "]]];
  ];

(* expected symbols per stage (common 01-03,05 + detector-specific 04,06) *)
syms01 = {"GeV","eV","keV","MeV","grams","cm","sec","Kel","me","Mpl",
          "km","mtr","nm","yr","month","day","kps","v0","ve","vesc",
          "alpha","\[Rho]DM","amukg","rhoAl","rhoTiN"};
syms02 = {"LIntpl1D","IntplArray","Inword"};
syms03 = {"\[Mu]\[Chi]t","ql","qr","q0","Jacobl","Jacobr","vmin","vmin0","Eq",
          "FDM","FmedH","FmedL","KKf","\[Eta]th","\[Eta]td","sigE","Eth"};
syms05 = {"\[Sigma]e","\[Sigma]N","TESsig","MKIDsig","TiNexp","Alexp"};

syms04TES = {"mp","AAl","mAl","a0Al","vLA","eps1AlIntpl","eps2AlIntpl",
             "ImepsAlIntpl","ImepsAlf","eps1Alf","eps2Alf","DwAl","Dw",
             "AveOmg","AveInvOmg","FnwAl","Fnw10GeV1","Fnw10GeV5"};
syms04MKID = {"kFTiN","vFTiN","wpTiN","GammaTiN","ImepsLTiN"};

syms06TES = {"energySum","midpointSum","BoxPDF","kerRAll","kerRAlr","RangeAlGen",
             "IntkerRAll","IntkerRAlr","KerRAll","KerRAlr",
             "CRTESLeft","CRTESRight","CRTES",
             "CRintTESLeft","CRintTESRight","CRintTES","CRintTESWeighted",
             "CRvmTES","CRvmTESLeft","CRvmTESRight"};
syms06MKID = {"energySum","midpointSum","BoxPDF","kerRTiNl","kerRTiNr","RangeTiNGen",
              "IntkerRTiNl","IntkerRTiNr","KerRTiNl","KerRTiNr",
              "CRTiNLeft","CRTiNRight","CRTiN",
              "CRintTiNLeft","CRintTiNRight","CRintTiN","CRintTiNWeighted",
              "CRvmTiN","CRvmTiNLeft","CRvmTiNRight"};

syms04 = If[det === "TES", syms04TES, syms04MKID];
syms06 = If[det === "TES", syms06TES, syms06MKID];

checkDefined["01_setup",             syms01];
checkDefined["02_functions_math",    syms02];
checkDefined["03_functions_response", syms03];
checkDefined["04_material",          syms04];
checkDefined["05_parameters",        syms05];
checkDefined["06_response_defs",     syms06];

(* ------------------------------------------------------ (C) numerical smoke *)
(* finiteRealQ: a genuine, finite, real number (machine or exact).            *)
finiteRealQ[x_] := NumericQ[x] && Element[x, Reals] && (Abs[N[x]] < Infinity);

smoke[label_String, expr_] :=
  Module[{val},
    val = Quiet[Check[expr, $Failed]];
    record["smoke/" <> label,
      If[finiteRealQ[val], "PASS", "FAIL"],
      If[finiteRealQ[val], ToString[N[val]], "got: " <> ToString[Short[val]]]];
  ];

(* common kinematics / form factor / halo *)
smoke["mu_chi_t",  \[Mu]\[Chi]t[1 GeV, me]];
smoke["FDM_n1",    FDM[100 keV][1]];
smoke["FDM_n0",    FDM[100 keV][0]];
smoke["vmin",      vmin[10 eV][1 GeV][100 keV]];
smoke["KKf",       KKf[v0, vesc]];
smoke["eta_th",    \[Eta]th[1 GeV][200 kps][v0, ve, vesc]];
smoke["energySum", energySum[(#^2) &, 1, 10, 0.5]];
smoke["midpointSum", midpointSum[(#) &, 0, 1, 10]];

If[det === "TES",
  (* material smoke (data-driven) *)
  smoke["ImepsAlf",  ImepsAlf[10., 1000.]];
  smoke["eps1Alf",   eps1Alf[10., 1000.]];
  smoke["Dw",        Dw[0.01]];
  smoke["AveOmg",    AveOmg];
  (* one full rate kernel evaluation (cheap point) *)
  smoke["kerRAll",   kerRAll[10 eV][1 GeV][1][220 kps]];
,
  (* MKID: analytic Lindhard dielectric *)
  smoke["ImepsLTiN", ImepsLTiN[1.0, 1000.]];
  smoke["wpTiN",     wpTiN];
  smoke["kerRTiNl",  kerRTiNl[10 eV][1 GeV][1][220 kps]];
];

(* --------------------------------------------------------------- reporting *)
nPass = Count[results, _?(#status === "PASS" &)];
nFail = Count[results, _?(#status === "FAIL" &)];
total = Length[results];

summary = StringJoin[
  "==== test_definitions report :: ", det, " ====\n",
  "loaded : ", detSrc, "\n",
  "passed : ", ToString[nPass], " / ", ToString[total], "\n",
  "failed : ", ToString[nFail], "\n",
  "----------------------------------------\n",
  StringRiffle[
    (StringJoin["[", #status, "] ", #name,
       If[#detail === "", "", "  -- " <> #detail]] & /@ results),
    "\n"],
  "\n"];

Export[FileNameJoin[{outDir, "report.txt"}], summary, "Text"];
Export[FileNameJoin[{outDir, "report.m"}],
  <|"detector" -> det, "passed" -> nPass, "failed" -> nFail,
    "total" -> total, "results" -> results|>];

Print["----------------------------------------"];
Print[det, ": ", nPass, "/", total, " passed, ", nFail, " failed"];
Print["report: ", FileNameJoin[{outDir, "report.txt"}]];

Exit[If[nFail == 0, 0, 1]];
