(* ::Package:: *)

(* ========================================================================== *)
(*  MKID - Natural-units eta(v_min) sampled on a response-matrix grid         *)
(*  TiN counterpart of TES/12_eta.wl.                                          *)
(*                                                                            *)
(*  The Python analysis contracts the response matrix M (natural units) with  *)
(*  a velocity distribution eta(v_min). For that product to be a real         *)
(*  expected event count, eta MUST be in the SAME natural units as M -- the   *)
(*  legacy data/Eta_data/*.csv files are in physical units (~cm^-1) and so    *)
(*  are NOT consistent with the new matrices. This stage produces the         *)
(*  consistent, natural-units eta:                                            *)
(*                                                                            *)
(*    eta_j = etaSHM[dmMass][ v_mid(j) ]      (natural units, GeV-based)       *)
(*                                                                            *)
(*  evaluated on the SAME v_min interval mid-points that label the matrix     *)
(*  columns (read from the matrix folder's vmin.csv). The result is written   *)
(*  next to matrix.csv so the Python loader can pick it up with no alignment   *)
(*  and no unit conversion.                                                    *)
(*                                                                            *)
(*  Physics: etaSHM = \[Eta]th from 03_functions_response.wl, the standard     *)
(*  halo (SHM) speed integral, which already carries the (rhoDM sigmae/mchi)   *)
(*  prefactor. CRTiN (06) carries NO rho/sigma/mass, so M @ eta does not        *)
(*  double-count. The SHM parameters v0, ve, vesc come from 01_setup.wl.       *)
(*                                                                            *)
(*  Scope: Halo (SHM) only. Disk/Bound are not defined in the current .wl      *)
(*  pipeline (03 has \[Eta]th == \[Eta]td and no Bound); their legacy CSVs are *)
(*  from an older pipeline and are intentionally not reproduced here.          *)
(*                                                                            *)
(*  Usage:                                                                     *)
(*      wolframscript -file 12_eta.wl <name|ALL> [model]                       *)
(*                                                                            *)
(*    <name>   substring of a matrix folder under output/MKID/response_matrix/ *)
(*             (e.g. TiN_q0M1_R5), or ALL for every folder found              *)
(*    [model]  velocity-distribution model; only "Halo" is implemented         *)
(*             (default Halo)                                                  *)
(*                                                                            *)
(*  e.g.  wolframscript -file 12_eta.wl ALL                                    *)
(*        wolframscript -file 12_eta.wl TiN_q0M1_R5 Halo                       *)
(*                                                                            *)
(*  Output (in each matched output/MKID/response_matrix/M*/<name>/ folder):    *)
(*    eta_<model>.csv   one column, length = #v_min columns of matrix.csv      *)
(* ========================================================================== *)

(* ============================ Arguments ============================ *)

commandLineArgs = Rest[$ScriptCommandLine];
nameArg  = If[Length[commandLineArgs] >= 1, commandLineArgs[[1]], "ALL"];
modelArg = If[Length[commandLineArgs] >= 2, Capitalize[ToLowerCase[commandLineArgs[[2]]]], "Halo"];

If[modelArg =!= "Halo",
  Print["error: only model 'Halo' is implemented (got ", modelArg, ")."];
  Print["       Disk/Bound are not defined in the current .wl pipeline."];
  Exit[1]];


(* ============================ Setup ============================ *)

fileDir   = DirectoryName[$InputFileName];
mathDir   = ParentDirectory[ParentDirectory[fileDir]];
matrixDir = FileNameJoin[{mathDir, "output", "MKID", "response_matrix"}];

(* Units + SHM params (01), eta definition + KKf (03, which Gets 02),
   reference cross sections sigmae (05, which Gets 04). *)
Get[FileNameJoin[{fileDir, "01_setup.wl"}]];
Get[FileNameJoin[{fileDir, "03_functions_response.wl"}]];
Get[FileNameJoin[{fileDir, "05_parameters.wl"}]];


(* ============================ Model / mass maps ============================ *)

(* Mass tag (M1/M2/M3) -> DM mass; matches constants.DM_MASS on the Python side
   and the "dmMass" stored per .wdx by 08. *)
dmMassOf[tag_] := Switch[tag,
  "M1", 10 MeV, "M2", 100 MeV, "M3", 1000 MeV, _, $Failed];

extractMass[name_] := Module[{m},
  m = StringCases[name, "M" ~~ d : DigitCharacter .. :> "M" <> d];
  If[m === {}, $Failed, First[m]]];

(* Natural-units halo speed integral evaluated at a v_min given in km/s.
   etaSHM[md][vKmS] = \[Eta]th[md][vKmS kps][v0, ve, vesc]. *)
etaSHM[md_][vKmS_] := \[Eta]th[md][vKmS kps][v0, ve, vesc];


(* ============================ Process one folder ============================ *)

processFolder[dir_] := Module[
  {vminFile, vminRows, vmid, tag, md, etaVals, outFile},

  vminFile = FileNameJoin[{dir, "vmin.csv"}];
  If[!FileExistsQ[vminFile], Print["skip (no vmin.csv): ", dir]; Return[]];

  vminRows = Import[vminFile, "CSV"];          (* rows of {v_low, v_high, v_mid} [km/s] *)
  vmid     = vminRows[[All, 3]];

  tag = extractMass[FileNameTake[dir]];
  md  = If[tag === $Failed, $Failed, dmMassOf[tag]];
  If[md === $Failed,
    Print["skip (no usable mass tag): ", dir]; Return[]];

  etaVals = N[etaSHM[md][#]] & /@ vmid;        (* natural units, one per v_min column *)

  outFile = FileNameJoin[{dir, "eta_" <> modelArg <> ".csv"}];
  Export[outFile, List /@ etaVals];            (* single column *)

  Print["Saved: ", FileNameTake[dir], "/eta_", modelArg, ".csv  (",
    Length[etaVals], " v_min points, ", modelArg, ", ", tag, ")"];
];


(* ============================ Dispatch ============================ *)

allFolders = DirectoryName /@ FileNames["vmin.csv", matrixDir, Infinity];

targets = If[ToUpperCase[nameArg] === "ALL",
  allFolders,
  Select[allFolders, StringContainsQ[#, nameArg] &]];

If[targets === {},
  Print["no matrix folders match ", nameArg, " under ", matrixDir],
  processFolder /@ targets];
