"""Physical constants -- the single source of truth for the Python side.

Every value is taken verbatim from the Mathematica pipeline so the two stay in
sync; the source file/line is given in a comment next to each. The convention
is natural units with everything expressed in GeV (``GeV = 1e9``), exactly as
in ``Mathematica/src/TES/01_setup.wl``.

Nothing here is tuned; downstream code derives its scale factors from these
named quantities instead of using bare literals.
"""

from __future__ import annotations

# ----- Energy base (01_setup.wl: "Energy base") -----------------------------
GeV = 1e9
# The Mathematica outputs (matrix.csv, eta_*.csv) were all computed at this
# GeV. Quantities loaded from them are rescaled from GEV_NATIVE to the active
# ``GeV`` by their energy-dimension in data_loader (matrix ~ GeV^-1, eta ~
# GeV^+1), so the Python unit system may differ from Mathematica's without
# re-running it. With GeV == GEV_NATIVE the rescaling is the identity.
GEV_NATIVE = 1e9
eV = 1e-9 * GeV
keV = 1e-6 * GeV
MeV = 1e-3 * GeV

# ----- Other base units in GeV equivalents (01_setup.wl) ---------------------
GRAM = 5.62e23 * GeV
CM = 1.0 / (1.98e-14 * GeV)
SEC = 1.0 / (6.58e-25 * GeV)

# ----- Mass scale (01_setup.wl) ---------------------------------------------
uG = 1e-6 * GRAM
NG = 1e-9 * GRAM
KG = 1e3 * GRAM
ME = 0.5109989 * MeV

# ----- Length / time (01_setup.wl) ------------------------------------------
KM = 1e5 * CM
MIC = 1e-4 * CM                # micrometre (Mathematica: mic = \[Mu]m)
NM = 1e-7 * CM                 # nanometre
KPS = KM / SEC                 # km/s, the unit of the v_min grid
YR = 365 * 24 * 3600 * SEC
MONTH = YR / 12
DAY = 24 * 3600 * SEC

# ----- Particle physics / DM density (01_setup.wl) --------------------------
ALPHA = 1.0 / 137
RHO_DM = 0.4 * GeV / CM**3

# ----- Material density (01_setup.wl) ---------------------------------------
RHO_AL = 2.7 * GRAM / CM**3
RHO_TIN = 5.4 * GRAM / CM**3

# ----- Reference cross sections (05_parameters.wl) --------------------------
SIGMA_E = 1e-30 * CM**2
SIGMA_N = 1e-32 * CM**2

# ----- Energy resolution (05_parameters.wl) ---------------------------------
TES_SIG = (0.068 / 0.8) / 2.355
MKID_SIG = 0.3 / 2.355

# ----- Detector exposures = active mass x time (05_parameters.wl) -----------
AL_EXP = 8200 * uG * MONTH      # Al (TES) exposure
# TiN (MKID): 1e7 detectors of TiN, over one year. Design based on arXiv: 2404.10785
TIN_EXP = 1e7 * 0.42 * NG * YR  # MKID (TiN) exposure

# Map the analysis mass tag (M1/M2/M3) to the DM mass; also stored per-file in
# each .wdx as "dmMass", so this table is only a convenience/cross-check.
DM_MASS = {"1": 10 * MeV, "2": 100 * MeV, "3": 1000 * MeV}
