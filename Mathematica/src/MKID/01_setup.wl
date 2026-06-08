(* ========================================================================== *)
(*  MKID - Setup (Constants, Units, common Parameters)                       *)
(*  Natural units with all quantities expressed in GeV.                         *)
(* ========================================================================== *)

fileDir = DirectoryName[$InputFileName];
Direc = FileNameJoin[{ParentDirectory[ParentDirectory[fileDir]], "input", "MKID"}];
SetDirectory[Direc];
Print["Direc = ", Direc];


(* ---- Energy base ---- *)
GeV = 10^9;
eV  = 10^(-9) GeV;
keV = 10^(-6) GeV;
MeV = 10^(-3) GeV;


(* ---- Other base units (GeV equivalents) ---- *)
grams = 5.62 * 10^23 GeV;
cm    = 1/(1.98 * 10^(-14) GeV);
sec   = 1/(6.58 * 10^(-25) GeV);
Kel   = 8.62 * 10^(-14) GeV;


(* ---- Mass scale ---- *)
ng = 10^(-9) grams;
\[Mu]g = 10^(-6) grams;
mg = 10^(-3) grams;
kg = 10^3 grams;
me = 0.5109989 MeV;
Mpl = 1.22 * 10^19 GeV;


(* ---- Length, Time ---- *)
km = 10^5 cm;
mtr = 10^2 cm;
mic = \[Mu]m = 10^(-4) cm;
nm = 10^(-7) cm;
yr = 365 * 24 * 3600 sec;
month = yr / 12;
day = 24 * 3600 sec;


(* ---- DM halo velocities (astro-ph/9710077v1, astro-ph/0611671) ---- *)
kps = km / sec;
v0 = 220 kps;
vemax = 241.279 kps;
vemin = 212.448 kps;
ve = (vemax + vemin) / 2;
vesc = 544 kps;


(* ---- Particle physics constants and DM density ---- *)
alpha = 1 / 137;
\[Rho]DM = 0.4 GeV / cm^3;
\[Rho]b = 10^14 GeV / cm^3;
amukg = 1.660538782 * 10^(-27);


(* ---- Material densities (used by 04_material) ---- *)
rhoAl = 2.7 grams / cm^3;
rhoTiN = 5.4 grams / cm^3;
