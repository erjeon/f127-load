"""화물-화물 회합 기하. 기하 인자(r^2, sin(theta))를 나눠서 저장한다.
   pipi.py TPR XTC RESN BEGIN_ps OUTDIR"""
import sys, numpy as np, MDAnalysis as mda
tpr, xtc, resn, beg, out = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), sys.argv[5]
u = mda.Universe(tpr, xtc)
g = u.select_atoms(f"resname {resn}")
res = list(g.residues)
if len(res) < 2: sys.exit("화물이 2개 미만이라 회합 분석을 건너뛴다")
de = np.arange(0, 1.601, 0.02); ae = np.arange(0, 90.1, 5.0)
hd = np.zeros(len(de)-1); ha = np.zeros(len(ae)-1)
def normal(p):
    q = p - p.mean(0)
    return np.linalg.svd(q)[2][2]
for ts in u.trajectory:
    if ts.time < beg: continue
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
print(f"  peak separation {cd[np.argmax(gd)]:.2f} nm ; parallel enrichment {ga[0]/ga[-1]:.1f}x")
print("  ★ raw 히스토그램은 기하 인자를 담고 있으므로 정규화한 3열을 쓸 것")
