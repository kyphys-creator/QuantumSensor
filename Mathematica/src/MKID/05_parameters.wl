(* ========================================================================== *)
(*  MKID - Calculation Parameters (cross sections, resolutions, exposure)    *)
(* ========================================================================== *)

fileDir = DirectoryName[$InputFileName];
Direc = FileNameJoin[{ParentDirectory[ParentDirectory[fileDir]], "input", "MKID"}];
Get[FileNameJoin[{fileDir, "04_material.wl"}]];
SetDirectory[Direc];
Print["Direc = ", Direc];


(* ---- Reference cross sections ---- *)
\[Sigma]e = 10^(-30) cm^2;
\[Sigma]N = 10^(-32) cm^2;


(* ---- Energy resolutions ---- *)
TESsig = (0.068 / 0.8) / 2.355;
MKIDsig = 0.3 / 2.355;


(* ---- Detector exposures (active mass times time) ---- *)
TiNexp = 10^7 rhoTiN (1 \[Mu]m) (50 \[Mu]m) (22 nm) month;
Alexp = (8200 ng) month;
