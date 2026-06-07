import numpy as np

from .constants import GRAM, YR, CM, SEC
from .data_loader import load_eta, load_response_matrix, load_ratebin
from .backgrounds import prepare_backgrounds, BKG_COLUMN
from .optimizer import run_optimize, run_optimize_qp
from .plotting import plot_eta_comparison


class NeutrinoAnalysis:
    """Neutrino flux optimization pipeline."""

    def __init__(self, background_scenario='c', eta='Halo', material='Al',
                 mass='1', q='0', cons1=1e-26, cons2=1e3):
        self.eta = eta
        self.material = material
        self.mass = mass
        self.q = q
        self.cons1 = cons1
        self.cons2 = cons2
        self.eps = 0

        # Notebook compatibility — expose unit factors as instance attrs.
        self.gram = GRAM
        self.cm = CM
        self.sec = SEC
        self.yr = YR

        self.fig1Solid = load_eta(eta, mass)
        self.CRmat = load_response_matrix(material, q, mass)
        self.Ratebin2 = load_ratebin(eta, material, q, mass)

        self.n = len(self.CRmat[0])
        self.m = len(self.Ratebin2)

        self.background_df = prepare_backgrounds()
        self.result = None
        self.set_background(background_scenario)

    def set_background(self, bkg_scenario='c'):
        self.background_scenario = bkg_scenario
        bkg_col = BKG_COLUMN.get(bkg_scenario)
        if bkg_col is None:
            raise ValueError(f"Invalid background scenario: {bkg_scenario}")

        self.ExtBkg = np.array(self.background_df[bkg_col]) / (1e3 * GRAM * YR)

        if self.background_scenario is 'none':
            self.M_matrix = self.cons1 * self.CRmat
            self.data_vector = self.cons2 * self.Ratebin2
            self.Bkg_vector = self.cons2 * self.ExtBkg
        else:
            self.M_matrix = self.cons1 * self.CRmat
            self.data_vector = self.cons2 * (self.Ratebin2 + self.ExtBkg)
            self.Bkg_vector = self.cons2 * self.ExtBkg

    def optimize(self, data=None, x0=None, display=True, solver='trust-constr',
                 fix=None, vertex_select=True):
        if data is None:
            data = self.data_vector
        if solver in ('osqp', 'qp', 'clarabel') or fix is not None:
            self.result = run_optimize_qp(
                self.M_matrix, data, self.Bkg_vector, self.n,
                eps=self.eps, fix=fix, vertex_select=vertex_select,
                verbose=display,
            )
        else:
            self.result = run_optimize(
                self.M_matrix, data, self.Bkg_vector, self.n,
                eps=self.eps, x0=x0, display=display,
            )
        return self.result

    def plot_eta_comparison(self, save=True):
        if self.result is None:
            raise RuntimeError("Please run optimization before plotting.")
        plot_eta_comparison(
            self.result.x, self.fig1Solid, self.n,
            self.cons1, self.cons2,
            self.eta, self.material, self.q, self.mass,
            self.background_scenario, save=save,
        )
