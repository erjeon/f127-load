"""Solute-solute association geometry, with the r^2 and sin(theta) weights divided out.
   pipi.py TPR XTC RESN BEGIN_ps OUTDIR"""
import sys, os, numpy as np, MDAnalysis as mda
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _universe import open_universe
tpr, xtc, resn, beg, out = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), sys.argv[5]
u = open_universe(tpr, xtc)
g = u.select_atoms(f"resname {resn}")
res = list(g.residues)
if len(res) < 2: sys.exit("fewer than two solute molecules, nothing to associate")
de = np.arange(0, 1.601, 0.02); ae = np.arange(0, 90.1, 5.0)
hd = np.zeros(len(de)-1); ha = np.zeros(len(ae)-1)
def normal(p):
    q = p - p.mean(0)
    return np.linalg.svd(q)[2][2]
nframes = 0
for ts in u.trajectory:
    if ts.time < beg: continue
    nframes += 1
    C = np.array([r.atoms.center_of_mass()/10.0 for r in res])
    Nv = np.array([normal(r.atoms.positions/10.0) for r in res])
    for i in range(len(res)):
        for j in range(i+1, len(res)):
            d = np.linalg.norm(C[i]-C[j])
            if d > de[-1]: continue
            a = np.degrees(np.arccos(min(1.0, abs(float(np.dot(Nv[i], Nv[j]))))))
            hd[np.searchsorted(de, d)-1] += 1
            ha[min(np.searchsorted(ae, a)-1, len(ha)-1)] += 1
cd = .5*(de[1:]+de[:-1]); ca = .5*(ae[1:]+ae[:-1])
gd = np.where(cd > 0, hd/cd**2, 0); ga = ha/np.sin(np.radians(ca))
np.savetxt(f"{out}/assoc_distance.dat", np.column_stack([cd, hd, gd/max(gd.max(),1e-9)]),
           header="r(nm)  raw_count  N(r)/r^2 (normalised)", fmt="%10.4f")
np.savetxt(f"{out}/assoc_angle.dat", np.column_stack([ca, ha, ga/max(ga.max(),1e-9)]),
           header="theta(deg)  raw_count  N(t)/sin(t) (normalised)", fmt="%10.4f")
# with no pair inside the outermost bin there is nothing to report
if hd.sum() == 0:
    print(f"  no {resn} pair came within {de[-1]:.1f} nm over these {nframes} frame(s), "
          "so there is no association to report")
else:
    enrich = ga[0] / ga[-1] if ga[-1] > 0 else float("inf")
    tail = f"{enrich:.1f}x" if np.isfinite(enrich) else "all pairs are face to face"
    print(f"  peak separation {cd[np.argmax(gd)]:.2f} nm, parallel enrichment {tail}")
    print("  the raw histogram carries the geometric weight. "
          "Use the normalised third column.")
