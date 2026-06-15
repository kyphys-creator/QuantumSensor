"""Verify colnorm weighting on the PRODUCTION (one-sided windowed) matrices,
without editing production. colnorm-in-x == uniform-in-z, so the proposed
one-line change is weight = D/D.max()  ->  weight = ones(n)."""
import numpy as np
from scipy.optimize import linprog
from quantum_sensor import DarkMatterQuantumAnalysis, RunConfig
from quantum_sensor.config import CONDITION_C
from quantum_sensor.model import condition
from quantum_sensor.optimizer import _build_qp, _column_scale

def colnorm_vertex(qp, D, mu):
    """Mirror optimizer._vertex_select but weight=ones in z (= colnorm in x)."""
    n = qp['n']; M_s = qp['M_active'] * D[None, :]
    A_ord_z = -(qp['A_ord'].toarray() * D[None, :])
    bounds = [(qp['eps']/d if d else 0.0, None) for d in D]
    res = linprog(np.ones(n), A_ub=A_ord_z, b_ub=np.zeros(n-1),
                  A_eq=M_s, b_eq=mu, bounds=bounds, method='highs-ds',
                  options={'presolve': True})
    return (D * res.x) if res.success else None

print(f"{'cfg':>16} {'uniform last':>13} {'colnorm last':>13}  {'unif err':>9} {'coln err':>9}")
for mass in ("1","2","3"):
    for nb in (5,10):
        a = DarkMatterQuantumAnalysis(RunConfig("Al","0",mass,nb,"Halo",background="none"))
        a.optimize()                              # production uniform vertex
        eta = a.eta
        m_cond,data,bkg,unscale = condition(a.m_phys,a.observed,a.background,c=CONDITION_C)
        qp = _build_qp(m_cond,data,bkg,a.n_vmin,0.0); D = _column_scale(qp['M_active'])
        xc = unscale(colnorm_vertex(qp, D, qp['y_active']))
        eu = np.linalg.norm(a.flux-eta)/np.linalg.norm(eta)
        ec = np.linalg.norm(xc-eta)/np.linalg.norm(eta)
        print(f"   Al M{mass} R{nb:<2} {a.flux[-1]/eta[-1]:13.2f} {xc[-1]/eta[-1]:13.2f}"
              f"  {eu:9.2e} {ec:9.2e}")
