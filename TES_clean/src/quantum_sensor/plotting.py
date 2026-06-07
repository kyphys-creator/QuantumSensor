from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from .constants import CM

RESULTS_DIR = Path(__file__).resolve().parents[2] / "results"

ENERGY_BIN_RANGES = {
    ('Al', '1'): (44, 535),
    ('Al', '2'): (13, 129),
    ('Al', '3'): (4, 128),
    ('TiN', '1'): (60, 259),
    ('TiN', '2'): (19, 249),
    ('TiN', '3'): (10, 248),
}


def save_optimized_flux(result_x, cons1, cons2, eta, material, q, mass):
    optimized_flux = result_x / cons2 * cons1 * CM
    out = RESULTS_DIR / f"Res5{eta}_{material}q{q}M{mass}.dat"
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(out, optimized_flux, delimiter=',')
    return optimized_flux, out


def plot_eta_comparison(result_x, fig1Solid, n, cons1, cons2,
                         eta, material, q, mass, background_scenario,
                         save=True):
    optimized_flux, _ = save_optimized_flux(result_x, cons1, cons2, eta, material, q, mass)

    lo, hi = ENERGY_BIN_RANGES[(material, mass)]
    optimized_energy_bins = np.linspace(lo, hi, n)

    x = np.linspace(1, 800, 800)
    plt.figure(figsize=(8, 6))
    plt.plot(x, fig1Solid, label=r'Expected', color='red')
    plt.hlines(optimized_flux[:-1], optimized_energy_bins[:-1], optimized_energy_bins[1:], linewidth=1.5)
    plt.xscale('log')
    plt.xlim(1, 800)
    plt.xlabel(r"$v_{min}$ [km/s]", fontsize=30)
    plt.ylabel(r"$\tilde{\eta}$ [cm$^{-1}$]", fontsize=30)
    plt.legend(fontsize=20)
    plt.grid(True, which="both", ls="--")

    if save:
        out_dir = RESULTS_DIR / f"scenario_bkg_{background_scenario}"
        out_dir.mkdir(parents=True, exist_ok=True)
        filename = out_dir / f"flux5_{background_scenario}_{material}M{mass}q{q}.pdf"
        plt.savefig(filename)
        print(f"Plot saved as {filename}")
