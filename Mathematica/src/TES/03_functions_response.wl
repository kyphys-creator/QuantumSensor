(* ========================================================================== *)
(*  TES - Domain / FDM / eta / Resolution functions                         *)
(* ========================================================================== *)

fileDir = DirectoryName[$InputFileName];
Direc = FileNameJoin[{ParentDirectory[ParentDirectory[fileDir]], "input", "TES"}];
Get[FileNameJoin[{fileDir, "02_functions_math.wl"}]];
SetDirectory[Direc];
Print["Direc = ", Direc];


(* ============================ Domain (kinematics) ============================ *)

(* Reduced mass of DM-target system *)
\[Mu]\[Chi]t[m\[Chi]_, mt_] := mt m\[Chi] / (mt + m\[Chi]);

(* Two kinematic branches of momentum transfer q at recoil energy Ee *)
ql[Ee_][m\[Chi]_][vmin_] := m\[Chi] vmin - Sqrt[m\[Chi]^2 vmin^2 - 2 m\[Chi] Ee];
qr[Ee_][m\[Chi]_][vmin_] := m\[Chi] vmin + Sqrt[m\[Chi]^2 vmin^2 - 2 m\[Chi] Ee];

(* Minimum q for given mass md (heavy mediator limit) *)
q0[Ee_][md_] := Sqrt[2 md Ee];

(* Jacobians dE/dq for left and right branches *)
Jacobl[Ee_][m\[Chi]_][vmin_] := -(m\[Chi] - m\[Chi]^2 vmin / Sqrt[m\[Chi]^2 vmin^2 - 2 m\[Chi] Ee]);
Jacobr[Ee_][m\[Chi]_][vmin_] :=  m\[Chi] + m\[Chi]^2 vmin / Sqrt[m\[Chi]^2 vmin^2 - 2 m\[Chi] Ee];

(* Minimum DM velocity *)
vmin[Ee_][m\[Chi]_][q_] := q / (2 m\[Chi]) + Ee / q;
vmin0[Ee_][md_] := Sqrt[2 Ee / md];

(* Energy-momentum relation for elastic scatter *)
Eq[q_][md_][v_] := q v - q^2 / (2 md);


(* ============================ FDM (DM form factor) ============================ *)

FDM[q_][n_] := (alpha me / q)^n;
FmedH[md_][q_] := ((md 240 kps)^2 + (100 md)^2) / (q^2 + (100 md)^2);
FmedL[md_][q_] := ((md 240 kps)^2 + 0) / (q^2 + 0);


(* ============================ eta (halo speed integrals) ============================ *)

KKf[v0_, vesc_] := v0^3 (-2. Exp[-vesc^2/v0^2] Pi (vesc/v0) + Pi^(3/2) Erf[vesc/v0]);

\[Eta]th[m\[Chi]_][vm_][v0_, ve_, vesc_] :=
  (\[Rho]DM \[Sigma]e / m\[Chi]) (v0^2 Pi / (2 ve KKf[v0, vesc])) Piecewise[{
    {-4 Exp[-vesc^2/v0^2] ve + Sqrt[Pi] v0 (Erf[(vm + ve)/v0] - Erf[(vm - ve)/v0]),
     vm < vesc - ve},
    {-2 Exp[-vesc^2/v0^2] (ve + vesc - vm) + Sqrt[Pi] v0 (Erf[vesc/v0] - Erf[(vm - ve)/v0]),
     vm < vesc + ve && vesc - ve < vm},
    {0, vm > vesc + ve}
  }];

\[Eta]td[m\[Chi]_][vm_][v0_, ve_, vesc_] :=
  (\[Rho]DM \[Sigma]e / m\[Chi]) (v0^2 Pi / (2 ve KKf[v0, vesc])) Piecewise[{
    {-4 Exp[-vesc^2/v0^2] ve + Sqrt[Pi] v0 (Erf[(vm + ve)/v0] - Erf[(vm - ve)/v0]),
     vm < vesc - ve},
    {-2 Exp[-vesc^2/v0^2] (ve + vesc - vm) + Sqrt[Pi] v0 (Erf[vesc/v0] - Erf[(vm - ve)/v0]),
     vm < vesc + ve && vesc - ve < vm},
    {0, vm > vesc + ve}
  }];


(* ============================ Resolution ============================ *)

sigE[Ep_, sig_] := sig Ep;
Eth[md_, vmax_] := md vmax^2 / 2;
