(* ========================================================================== *)
(*  MKID - Material (TiN, analytic Lindhard dielectric function)                *)
(* ========================================================================== *)

fileDir = DirectoryName[$InputFileName];
Direc = FileNameJoin[{ParentDirectory[ParentDirectory[fileDir]], "input", "MKID"}];
Get[FileNameJoin[{fileDir, "03_functions_response.wl"}]];
SetDirectory[Direc];
Print["Direc = ", Direc];


(* ============================ TiN constants ============================ *)

rhoTiN = 5.4 grams / cm^3;
kFTiN = Sqrt[2 me (3.94 eV)];
vFTiN = kFTiN / me;
wpTiN = Sqrt[(4 Pi alpha / me) (2 me (3.94 eV))^(3/2) / (3 Pi^2)];
GammaTiN = 0.1 (3.94 eV);


(* ============================ TiN dielectric (Lindhard) ============================ *)

ImepsLTiN[w_, q_] := Im[-((1 + (3 wpTiN^2 / (q^2 vFTiN^2)) *
  (1/2 + (kFTiN / (4 q)) *
    (1 - ((q / (2 kFTiN)) - (w + I GammaTiN) / (q vFTiN))^2) *
    Log[((q / (2 kFTiN)) - (w + I GammaTiN) / (q vFTiN) + 1) /
        ((q / (2 kFTiN)) - (w + I GammaTiN) / (q vFTiN) - 1)] +
    (kFTiN / (4 q)) *
    (1 - ((q / (2 kFTiN)) + (w + I GammaTiN) / (q vFTiN))^2) *
    Log[((q / (2 kFTiN)) + (w + I GammaTiN) / (q vFTiN) + 1) /
        ((q / (2 kFTiN)) + (w + I GammaTiN) / (q vFTiN) - 1)]))^(-1))];
