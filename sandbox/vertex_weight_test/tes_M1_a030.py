import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from weight_test import (central_area_cut, trap_weights, eta_halo_at,
                         weighted_min_flux_lp, CM_at, GEV_NATIVE, CSV_DIR, TARGET)

ALPHA=0.30; mass="1"; GEV=1e9/ \
    abs(eta_halo_at(GEV_NATIVE,"1",np.genfromtxt(CSV_DIR/"Al_q0M1_R5.csv",delimiter=",",names=True)["vmin"])).max()
d=np.genfromtxt(CSV_DIR/"Al_q0M1_R5.csv",delimiter=",",names=True); v=d["vmin"]
R=np.array([d[c] for c in d.dtype.names if c!="vmin"]); w=trap_weights(v)
Rc=np.array([central_area_cut(v,R[i],ALPHA) for i in range(len(R))])
M0=(Rc*w[None,:]); keep=M0.any(axis=0); M0k,vk=M0[:,keep],v[keep]
M=M0k/(GEV/GEV_NATIVE); eta=eta_halo_at(GEV,mass,vk)
mp=(TARGET/float((M@eta).max()))*M; data=mp@eta
cn=np.linalg.norm(mp,axis=0)
C=CM_at(GEV)
xu=weighted_min_flux_lp(mp,data,np.ones(mp.shape[1]))
xc=weighted_min_flux_lp(mp,data,cn/cn.max())
print(f"window [{vk.min():.0f},{vk.max():.0f}] km/s, {len(vk)} cols")
print(f"uniform: last x/eta={xu[-1]/eta[-1]:.3f}  err={np.linalg.norm(xu-eta)/np.linalg.norm(eta):.2e}")
print(f"colnorm: last x/eta={xc[-1]/eta[-1]:.3f}  err={np.linalg.norm(xc-eta)/np.linalg.norm(eta):.2e}")
vg=np.logspace(0,np.log10(800),400)
fig,ax=plt.subplots(figsize=(8,5.5))
ax.plot(vg, eta_halo_at(GEV,mass,vg)*C, color="red", lw=2.5, label="input eta (SHM)")
ax.step(vk, xu*C, where="mid", color="C1", lw=1.8, label=f"uniform (last {xu[-1]/eta[-1]:.2f})")
ax.step(vk, xc*C, where="mid", color="C0", lw=1.8, label=f"colnorm (last {xc[-1]/eta[-1]:.2f})")
ax.axvspan(vk.min(),vk.max(),color="C0",alpha=0.06,label="matrix window")
ax.set_xscale("log"); ax.set_xlim(1,800); ax.set_ylim(bottom=0)
ax.set_xlabel(r"$v_{min}$ [km/s]"); ax.set_ylabel(r"flux $\tilde{\eta}$ [cm$^{-1}$]")
ax.set_title("TES Al q0 M1 R5 SHM, alpha=0.3  (uniform vs colnorm)")
ax.grid(True,which="both",ls="--",alpha=0.35); ax.legend(fontsize=10)
fig.savefig("tes_M1_a030.pdf",bbox_inches="tight"); print("saved tes_M1_a030.pdf")
