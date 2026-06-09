"""Minimal end-to-end example.

Runs the forward+inverse analysis for Al, heavy mediator (q0), 5 energy bins,
Halo halo model, signal only -- for the three DM masses -- and saves a
flux-comparison figure for each.

    python examples/run_example.py
"""

from quantum_sensor import DarkMatterQuantumAnalysis, RunConfig

for mass in ("1", "2", "3"):
    cfg = RunConfig(material="Al", q="0", mass=mass, nbins=5,
                    eta="Halo", background="none")
    analysis = DarkMatterQuantumAnalysis(cfg)
    flux = analysis.optimize(solver="osqp")
    print(f"M{mass}: {analysis.n_ebins} energy bins x {analysis.n_vmin} v_min bins,"
          f" peak flux {flux.max():.3e}")
    analysis.save_flux()
    analysis.plot(save=True)
