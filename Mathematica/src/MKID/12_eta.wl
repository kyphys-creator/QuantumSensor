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
(*  Physics: eta = \[Eta]th from 03_functions_response.wl, the truncated-       *)
(*  Maxwellian speed integral, which already carries the (rhoDM sigmae/mchi)    *)
(*  prefactor. CRTiN (06) carries NO rho/sigma/mass, so M @ eta does not        *)
(*  double-count.                                                              *)
(*                                                                            *)
(*  Models:                                                                    *)
(*    Halo   standard halo (SHM): \[Eta]th with v0, ve, vesc (rhoDM)            *)
(*    Disk   pure dark disk: \[Eta]th with v0DD, veDD, vescDD (same rhoDM)      *)
(*    Bound  Earth-bound thermal DM: a DIFFERENT, dense (rho_b) isotropic        *)
(*           (ve=0) population; isotropic truncated Maxwellian with thermal      *)
(*           v0 = Sqrt[2 kT/mchi] (~2.16 km/s at 1 GeV) and vesc = 11.2 km/s.    *)
(*           Nonzero only for v_min < 11.2 km/s, so in practice only the         *)
(*           heaviest mass (M3) whose window reaches that low has any signal.    *)
(*                                                                            *)
(*  Usage:                                                                     *)
(*      wolframscript -file 12_eta.wl <name|ALL> [model]                       *)
(*                                                                            *)
(*    <name>   substring of a matrix folder under output/MKID/response_matrix/ *)
(*             (e.g. TiN_q0M1_R5), or ALL for every folder found              *)
(*    [model]  velocity model: Halo (default), Disk, or Bound                  *)
(*                                                                            *)
(*  e.g.  wolframscript -file 12_eta.wl ALL                                    *)
(*        wolframscript -file 12_eta.wl ALL Disk                               *)
(*        wolframscript -file 12_eta.wl ALL Bound                              *)
(*        wolframscript -file 12_eta.wl TiN_q0M1_R5 Halo                       *)
(*                                                                            *)
(*  Output (in each matched output/MKID/response_matrix/M*/<name>/ folder):    *)
(*    eta_<model>.csv   one column, length = #v_min columns of matrix.csv      *)
(* ========================================================================== *)

(* ============================ Arguments ============================ *)

commandLineArgs = Rest[$ScriptCommandLine];
nameArg  = If[Length[commandLineArgs] >= 1, commandLineArgs[[1]], "ALL"];
modelArg = If[Length[commandLineArgs] >= 2, Capitalize[ToLowerCase[commandLineArgs[[2]]]], "Halo"];


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

(* Velocity-distribution models. Halo and Disk share 03's truncated-Maxwellian
   \[Eta]th (same DM density rhoDM, only the {v0, ve, vesc} velocities differ):
     Halo = standard halo (SHM):       v0, ve, vesc       from 01_setup
     Disk = pure dark disk (colder):   v0DD, veDD, vescDD  (same DM density)

   Bound = Earth-bound thermalized DM: a DIFFERENT, much denser population
   (rho_b = 1e14 GeV/cm^3) that is isotropic (ve = 0). The ve -> 0 limit of
   \[Eta]th is an isotropic truncated Maxwellian, defined directly here to avoid
   the 1/ve in \[Eta]th. Its dispersion is thermal, v0 = Sqrt[2 kT / mchi]
   (~2.16 km/s at 1 GeV for kT = 300 K), and vesc = 11.2 km/s. *)
TEB    = 300 Kel;          (* Earth temperature; sets the thermal v0 below *)
vescEB = 11.2 kps;
v0EB[md_] := Sqrt[2 TEB / md];
\[Eta]bound[md_][vm_] := With[{v0 = v0EB[md]},
  If[vm >= vescEB, 0,
    (\[Rho]b \[Sigma]e / md) (2 Pi v0^2 / KKf[v0, vescEB]) *
      (Exp[-vm^2/v0^2] - Exp[-vescEB^2/v0^2])]];

If[!MemberQ[{"Halo", "Disk", "Bound"}, modelArg],
  Print["error: model must be 'Halo', 'Disk', or 'Bound' (got ", modelArg, ")."];
  Exit[1]];

(* Natural-units speed integral at a v_min given in km/s, for the chosen model. *)
etaModel[md_][vKmS_] := Switch[modelArg,
  "Halo",  \[Eta]th[md][vKmS kps][v0, ve, vesc],
  "Disk",  \[Eta]th[md][vKmS kps][v0DD, veDD, vescDD],
  "Bound", \[Eta]bound[md][vKmS kps]];


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

  etaVals = N[etaModel[md][#]] & /@ vmid;       (* natural units, one per v_min column *)

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
