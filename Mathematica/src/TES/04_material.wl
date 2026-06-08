(* ========================================================================== *)
(*  TES - Material (Al, Mermin model dielectric function)                       *)
(* ========================================================================== *)

fileDir = DirectoryName[$InputFileName];
Direc = FileNameJoin[{ParentDirectory[ParentDirectory[fileDir]], "input", "TES"}];
Get[FileNameJoin[{fileDir, "03_functions_response.wl"}]];
SetDirectory[Direc];
Print["Direc = ", Direc];


(* ============================ Al constants ============================ *)

rhoAl = 2.7 grams / cm^3;
mp = 0.94 GeV;
AAl = 27;
mAl = mp AAl;
a0Al = 4.046 * 10^(-8) cm;
vLA = 6.2 kps;


(* ============================ Al dielectric (Mermin model) ============================ *)

datAl = Import["Al_mermin.dat"];
datLAl = Length[datAl];

eps1Aldat  = ConstantArray[0, datLAl - 1];
eps2Aldat  = ConstantArray[0, datLAl - 1];
ImepsAldat = ConstantArray[0, datLAl - 1];

For[i = 2, i <= datLAl, i++,
  eps1Aldat[[i - 1]]  = {{datAl[[i, 1]], datAl[[i, 2]]}, datAl[[i, 3]]};
  eps2Aldat[[i - 1]]  = {{datAl[[i, 1]], datAl[[i, 2]]}, datAl[[i, 4]]};
  ImepsAldat[[i - 1]] = {{datAl[[i, 1]], datAl[[i, 2]]},
                         datAl[[i, 4]] / (datAl[[i, 4]]^2 + datAl[[i, 3]]^2)};
];

eps1AlIntpl  = Interpolation[eps1Aldat,  InterpolationOrder -> 1];
eps2AlIntpl  = Interpolation[eps2Aldat,  InterpolationOrder -> 1];
ImepsAlIntpl = Interpolation[ImepsAldat, InterpolationOrder -> 1];

ImepsAlf[w_, q_] := If[(0 <= w <= 99.3) && (0 <= q <= 37289.5), ImepsAlIntpl[w, q], 0];
eps1Alf[w_, q_]  := If[(0 <= w <= 99.3) && (0 <= q <= 37289.5), eps1AlIntpl[w, q],  0];
eps2Alf[w_, q_]  := If[(0 <= w <= 99.3) && (0 <= q <= 37289.5), eps2AlIntpl[w, q],  0];


(* ============================ Al phonon density of states ============================ *)

AlRawDos = Transpose[Delete[Import[FileNameJoin[{Direc, "Al_pDoS.dat"}]], 1]];
Print[{AlRawDos[[1, 1]], AlRawDos[[1, -1]]}];
DwAl = IntplArray[AlRawDos];
wDw = AlRawDos[[1, -1]];

Dw[w_?NumericQ] := Piecewise[{{DwAl[[1]][w], w <= wDw}, {0, w > wDw}}];

AveOmg    = NIntegrate[w     DwAl[[1]][w], {w, 0., wDw}];
AveInvOmg = NIntegrate[DwAl[[1]][w] / w,   {w, 0., wDw}];


(* ============================ Al form factor table ============================ *)

AlRawFn = Transpose[Delete[Import[FileNameJoin[{Direc, "Al_Fn.dat"}]], 1]];
Print[{AlRawFn[[1, 1]], AlRawFn[[1, -1]]}];
FnwAl = IntplArray[AlRawFn];

Fnw10GeV1[w_] := Piecewise[{{FnwAl[[1]][w], w <= 0.2}, {0, w > 0.2}}];
Fnw10GeV2[w_] := Piecewise[{{FnwAl[[2]][w], w <= 0.2}, {0, w > 0.2}}];
Fnw10GeV3[w_] := Piecewise[{{FnwAl[[3]][w], w <= 0.2}, {0, w > 0.2}}];
Fnw10GeV4[w_] := Piecewise[{{FnwAl[[4]][w], w <= 0.2}, {0, w > 0.2}}];
Fnw10GeV5[w_] := Piecewise[{{FnwAl[[5]][w], w <= 0.2}, {0, w > 0.2}}];
