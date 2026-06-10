(* ========================================================================== *)
(*  TES - Calculation Parameters (cross sections, resolutions, exposure)    *)
(* ========================================================================== *)

fileDir = DirectoryName[$InputFileName];
Direc = FileNameJoin[{ParentDirectory[ParentDirectory[fileDir]], "input", "TES"}];
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
(* TiN (MKID): 10^7 detectors, 0.42 ng each, over one year. Design based on arXiv:2404.10785. *)
TiNexp = 10^7 (0.42 ng) yr;
Alexp = (8200 \[Mu]g) month;
