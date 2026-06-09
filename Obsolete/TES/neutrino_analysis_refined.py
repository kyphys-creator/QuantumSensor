import numpy as np
from scipy.optimize import minimize, LinearConstraint
from scipy import integrate
from scipy.sparse import csc_matrix
import time
import matplotlib.pyplot as plt
import pandas as pd
from numba import jit, njit
from matplotlib.ticker import LogLocator

class NeutrinoAnalysis:
    """
    A class to encapsulate the entire neutrino flux optimization analysis.

    This class handles data loading, constants definition, optimization,
    parameter scanning, and plotting. It is designed to be modular, allowing
    different background scenarios to be tested easily.
    """
    def __init__(self, background_scenario='c', eta = 'Halo', material = 'Al', mass = '1', q= '0', cons1 = 1e-26, cons2 = 1e3):
        """
        Initializes the analysis pipeline.

        This constructor acts as the main entry point for the class, orchestrating
        the setup of the analysis environment.

        Args:
            background_scenario (str): A character indicating the background model
                                       to use. Accepted values are 'a', 'b', 'c',
                                       or 'flat'. Defaults to 'c'.
        """
        self.iterationTime = 1000  # Number of points for integration steps

        # --- Setup Sequence ---
        self._load_data(eta,material,mass,q)              # Load all necessary data from files.
        self._define_constants(cons1,cons2)       # Define physical and numerical constants.
        self._prepare_backgrounds()    # Pre-calculate different background models.
        
        # This will be populated by the main optimization result.
        self.result = None
        
        # Set the specified background model and prepare matrices for optimization.
        self.set_background(background_scenario)

    def _load_data(self, eta = 'Halo', material = 'Al', mass = '1', q = '0'):
        """
        Loads all required data from external CSV files.

        This includes theoretical flux data from "Danny's files", the detector
        response matrix (CRmat), and event rate data.
        """
        self.eta = eta
        self.material = material
        self.mass = mass
        self.q = q
        # Load theoretical flux data for comparison plots.
        self.fig1Solid = np.genfromtxt(f'Eta_data/Eta{self.eta}M{self.mass}_Ko.csv', delimiter=',')
        
        # Load experimental/simulation data.
        # CRmat is the detector response matrix, mapping true flux to observed events.


        #self.CRmat = np.genfromtxt(f'CRmat/CR5mat{self.material}q{self.q}M{self.mass}_Ko_Refined.csv', delimiter=',')# in natural unit GeV = 1.
        self.CRmat = np.genfromtxt(f'curlyRplotting/{self.material}_R5_q{self.q}M{self.mass}_Ko_thresholded_refined_exposure.csv', delimiter=',')# in natural unit GeV = 1.
        # Ratebin7 and Ratebin2 are likely event rates from different energy bins or time periods.
        #self.Ratebin2 = np.genfromtxt(f'Ratebin/Event5Eta{self.eta}Mat{self.material}q{self.q}M{self.mass}.csv', delimiter=',')# in natural unit GeV = 1.
        self.Ratebin2 = np.genfromtxt(f'Ratebin/Event5Eta{self.eta}Mat{self.material}q{self.q}M{self.mass}.csv', delimiter=',')# in natural unit GeV = 1.
        # RateDiff is the difference in rates, which is part of the background model.

    def _define_constants(self, cons1 = 1e-26, cons2 = 1e3):
        """
        Defines physical and numerical constants used throughout the analysis.
        """
        # --- Dimensions ---
        # n: Number of parameters to optimize (dimension of the flux vector).
        self.n = len(self.CRmat[0])
        # m: Number of data points/energy bins.
        self.m = len(self.Ratebin2)
        
        # --- Physical Constants (in specific unit system) ---
        GeV = 1e9
        self.gram = 5.62e23 * GeV    # Conversion factor for grams.
        self.cm = (1.98e-14 * GeV)**(-1)  # Conversion factor for centimeters.
        self.sec = (6.58e-25 * GeV)**(-1)   # Conversion factor for seconds.
        self.yr = np.pi * 1e7 * self.sec   # Conversion factor for years.

        self.cons1=cons1
        self.cons2=cons2

        # --- Numerical Constants ---
        # A small value to prevent division by zero or taking the log of zero.
        self.eps = 0
        
        # --- Detector and Binning ---
        # Defines the energy bin edges in eV for background calculations.
        self.bins = np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1])
        # Mass of the detector in kg.
        self.Mdetector = 1.0 #kg

    def _compute_bi(self, A, B, C):
        """
        Helper function to compute background event rates per bin.

        This calculates the integral of a background model of the form:
        R(E) = A * exp(-E/B) + C
        across each energy bin.

        Args:
            A (float): Normalization of the exponential component.
            B (float): 'Decay' constant of the exponential component.
            C (float): Constant or 'flat' component.

        Returns:
            np.array: An array of background counts for each bin.
        """
        bi = []
        for i in range(len(self.bins) - 1):
            Ei, Ei1 = self.bins[i], self.bins[i+1]
            # Analytical integral of the model over the bin [Ei, Ei1].
            integral = -A * B * (np.exp(-Ei1/B) - np.exp(-Ei/B)) + C * (Ei1 - Ei)
            bi.append(self.Mdetector * integral)
        return np.array(bi)

    def _prepare_backgrounds(self):
        """
        Computes and stores different background scenarios in a pandas DataFrame.

        This method pre-calculates the background rates for several plausible
        scenarios (a, b, c, flat) for easy access later.
        """
        eV = 1
        keV = 1000
        days_in_year = 365
        
        # Define parameters for each scenario.
        # Scenario (a): Zero background.
        ba = self._compute_bi(A=0 * eV / keV * days_in_year, B=10, C=0 * eV / keV * days_in_year)
        # Scenario (b): High background.
        bb = self._compute_bi(A=460 * eV / keV * days_in_year, B=10, C=100 * eV / keV * days_in_year)
        # Scenario (b)x2: High background.
        b2 = self._compute_bi(A=920 * eV / keV * days_in_year, B=10, C=200 * eV / keV * days_in_year)
        # Scenario (c): Intermediate background.
        bc = self._compute_bi(A=50 * eV / keV * days_in_year, B=10, C=20 * eV / keV * days_in_year)
        # Scenario (flat): Constant background only.
        bflat = self._compute_bi(A=0 * eV / keV * days_in_year, B=10, C=100 * eV / keV * days_in_year)

        # Store the computed backgrounds in a DataFrame.
        self.background_df = pd.DataFrame({
            "Bin Start [eV]": self.bins[:-1],
            "Bin End [eV]": self.bins[1:],
            "b_i (a)": ba,
            "b_i (b)": bb,
            "b_i (c)": bc,
            "b_i (flat)": bflat,
            "b_i (b2)": b2,
        })
        
    def set_background(self, bkg_scenario='c'):
        """
        Sets the external background for the analysis and pre-computes matrices.

        This is a key method that prepares the final vectors and matrices needed
        for the chi-squared calculation based on the chosen background scenario.
        Pre-computing these avoids redundant calculations inside the optimizer loop.
        
        Args:
            bkg_scenario (str): The background scenario to apply ('a', 'b', 'c', 'flat', 'none').
        """
        self.background_scenario = bkg_scenario
        bkg_map = {
            'a': 'b_i (a)', 'b': 'b_i (b)',
            'c': 'b_i (c)', 'flat': 'b_i (flat)', 'none': 'b_i (a)',
            'b2': 'b_i (b2)'
        }
        bkg_col = bkg_map.get(bkg_scenario)
        if bkg_col is None:
            raise ValueError(f"Invalid background scenario: {bkg_scenario}")

        # Select the external background vector and scale it by the appropriate units.
        self.ExtBkg = np.array(self.background_df[bkg_col]) / (1e3 * self.gram * self.yr) #1/(kg yrs) -> GeV = 1 natural unit
        
        # This multiplier scales the entire equation for numerical stability.
        # --- Pre-compute core components for optimization ---
        # M_matrix: The full response matrix, scaled.
        if self.background_scenario is 'none':
            self.M_matrix = self.cons1 * self.CRmat
            self.data_vector = self.cons2 * self.Ratebin2
            self.Bkg_vector = self.cons2 * (self.ExtBkg)
        else:
            self.M_matrix = self.cons1 * self.CRmat
        # data_vector: The observed data vector, including the chosen background.
            self.data_vector = self.cons2 * (self.Ratebin2 + self.ExtBkg)
        # Bkg_vector: The total background model to be subtracted from the data.
            self.Bkg_vector = self.cons2 * (self.ExtBkg)

    @staticmethod
    @njit
    def _chiN(xvec, M, data, Bkg, eps):
        """
        Calculates the Neyman chi-squared value.

        This static method is decorated with @njit from the Numba library, which
        compiles it to fast native machine code for significant speedup.

        chi^2 = sum( (data - Bkg - model)^2 / data )

        Args:
            xvec (np.array): The vector of parameters being optimized (the flux).
            M (np.array): The scaled response matrix.
            data (np.array): The scaled data vector.
            Bkg (np.array): The scaled background vector.
            eps (float): A small number to prevent division by zero.

        Returns:
            float: The calculated chi-squared value.
        """
        # Calculate the expected number of events from the flux: model = M * xvec
        mod = M.dot(xvec)
        # Ensure the model prediction is positive to avoid issues in the denominator.
        #mod = np.maximum(mod, eps)

        mask = data > 0
        # Neyman's chi-squared formula.
        return np.sum(((data[mask] - Bkg[mask] - mod[mask]) ** 2) / data[mask])
# 67.8%
# 95.4%
    def optimize(self, data, x0=None, display = True):
        """
        Runs the main optimization to find the best-fit flux vector.

        This method sets up the constraints and calls the minimizer.

        Args:
            x0 (np.array, optional): Initial guess for the parameters. If None,
                                     a vector of ones is used.

        Returns:
            scipy.optimize.OptimizeResult: The result object from the minimizer.
        """
        if x0 is None:
            x0 = np.ones(self.n) # Default initial guess.
            
        # --- Define Constraints ---
        # 1. Ordering Constraint: Ensures the flux is non-increasing with energy.
        #    This means x[i] >= x[i+1] for all i.
        A = np.zeros((self.n - 1, self.n))
        for i in range(self.n - 1):
            A[i, i] = 1
            A[i, i + 1] = -1
        A_sparse = csc_matrix(A) # Use sparse matrix for efficiency.
        ordering_constraint = LinearConstraint(A_sparse, lb=0, ub=np.inf)

        # 2. Bounds: Ensures the flux is strictly positive.
        bounds = [(self.eps, None)] * self.n
        def hessian_chiN(xvec):
            safe_data = np.where(data > 0, data, 1)
            weighted = (2 * self.M_matrix.T / safe_data) @ self.M_matrix
            return weighted
        # --- Run Minimization ---
        if display:
            res = minimize(
                # Objective function to minimize (our chi-squared).
                lambda x: self._chiN(x, self.M_matrix, data, self.Bkg_vector, self.eps),
                x0,
                method='trust-constr',  # A trust-region method suitable for constraints.
                bounds=bounds,
                constraints=[ordering_constraint],
                options={'maxiter': 100000, 'xtol': 1e-9, 'gtol': 1e-9,'verbose': 3},
                jac = '3-point',
                #hess=hessian_chiN
            )
        else:
            res = minimize(
                # Objective function to minimize (our chi-squared).
                lambda x: self._chiN(x, self.M_matrix, data, self.Bkg_vector, self.eps),
                x0,
                method='trust-constr',  # A trust-region method suitable for constraints.
                bounds=bounds,
                constraints=[ordering_constraint],
                options={'maxiter': 100000, 'xtol': 1e-9, 'gtol': 1e-9},
                jac = '3-point',
                #hess=hessian_chiN
            )
        self.result = res # Store the result within the class instance.
        return res
        
    
    def plot_eta_comparison(self, save=True):
        """
        Plots the final optimized flux against the theoretical integrated flux data.

        Args:
            save (bool): If True, saves the plot to a PDF file.
        """
        if self.result is None:
            raise RuntimeError("Please run optimization before plotting.")
        
        # Get the theoretical curves.
        x=np.linspace(1, 800, 800)
        Phi=self.fig1Solid

        plt.figure(figsize=(8, 6))
        # Plot theoretical curves.
        plt.plot(x, Phi, label = r'Expected', color='red')
        # Plot the discrete flux points found by the optimization.
        if self.material == 'Al':
            if self.mass == '1':
                optimized_energy_bins = np.linspace(44, 535, self.n)
            elif self.mass == '2':
                optimized_energy_bins = np.linspace(13, 129, self.n)
            elif self.mass == '3':
                optimized_energy_bins = np.linspace(4, 128, self.n)
        elif self.material == 'TiN':
            if self.mass == '1':
                optimized_energy_bins = np.linspace(60, 259, self.n)
            elif self.mass == '2':
                optimized_energy_bins = np.linspace(19, 249, self.n)
            elif self.mass == '3':
                optimized_energy_bins = np.linspace(10, 248, self.n)
        # Scale the optimized parameters back to physical flux units.
        optimized_flux = self.result.x/self.cons2 * self.cons1 *self.cm
        np.savetxt(f'Res5{self.eta}_{self.material}q{self.q}M{self.mass}.dat', optimized_flux, delimiter=',')
        #plt.scatter(optimized_energy_bins, optimized_flux, label=r'$\chi^2$ minimization', zorder=5, s = 1)
        plt.hlines(optimized_flux[:-1], optimized_energy_bins[:-1], optimized_energy_bins[1:], linewidth=1.5)
        
        plt.xscale('log')
        #plt.ylim(0, 1e-30)
        plt.xlim(1, 800)
        plt.xlabel(r"$v_{min}$ [km/s]", fontsize = 30)
        plt.ylabel(r"$\tilde{\eta}$ [cm$^{-1}$]", fontsize = 30)
        #plt.title('Halo function')
        plt.legend(fontsize = 20)
        plt.grid(True, which="both", ls="--")
        
        if save:
            filename = f'scenario_bkg_{self.background_scenario}/flux5_{self.background_scenario}_{self.material}M{self.mass}q{self.q}.pdf'
            plt.savefig(filename)
            print(f"Plot saved as {filename}")
        #plt.show()