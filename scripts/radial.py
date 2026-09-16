"""Radial density, uptake and solute hydration.  radial.py TPR XTC RESN BEGIN_ps OUTDIR"""
import sys, os, numpy as np, MDAnalysis as mda
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _universe import open_universe
tpr, xtc, resn, beg, out = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), sys.argv[5]
u = open_universe(tpr, xtc)
core   = u.select_atoms("resname PROXS PROXR PROX")
corona = u.select_atoms("resname ETHOX ETHO")
guest  = u.select_atoms(f"resname {resn}")
water  = u.select_atoms("resname TIP3 SOL and name OH2 OW O")
mic = core + corona
assert len(core) and len(guest), "the core or solute selection is empty"
edges = np.linspace(0, 7.0, 29); cent = .5*(edges[1:]+edges[:-1])
V = (4/3)*np.pi*(edges[1:]**3-edges[:-1]**3)
acc = {k: np.zeros(len(cent)) for k in ('core','corona','guest','water')}
R90, enc, nwat, dg, n = [], [], [], [], 0
for ts in u.trajectory:
    if ts.time < beg: continue
    com = mic.center_of_mass()/10.0
    rad = lambda g: np.linalg.norm(g.positions/10.0 - com, axis=1)
    rc = rad(core)
    acc['core']   += np.histogram(rc, bins=edges)[0]/V
    acc['corona'] += np.histogram(rad(corona), bins=edges)[0]/V
    acc['guest']  += np.histogram(rad(guest),  bins=edges)[0]/V
    acc['water']  += np.histogram(rad(water),  bins=edges)[0]/V
    r90 = np.percentile(rc, 90); R90.append(r90)
    dres = np.array([np.linalg.norm(r.atoms.center_of_mass()/10.0-com) for r in guest.residues])
    dg.extend(dres); enc.append((ts.time/1000.0, int((dres < r90).sum())))
    gp = guest.positions/10.0; wp = water.positions/10.0
    dd = wp[None,:,:]-gp[:,None,:]
    nwat.append((ts.time/1000.0, int((np.linalg.norm(dd,axis=2) < 0.35).sum())))
    n += 1
for k in acc: acc[k] /= n
R90 = float(np.mean(R90)); dg = np.array(dg)
np.savetxt(f"{out}/radial_density.dat",
           np.column_stack([cent, acc['core'], acc['corona'], acc['guest'], acc['water']]),
           header=f"r(nm) core corona {resn} water   ; frames={n}  R90={R90:.3f} nm", fmt="%10.4f")
np.savetxt(f"{out}/encapsulation.dat", np.array(enc), header=f"time(ns) n_inside_R90 (total {len(guest.residues)})", fmt="%10.3f")
np.savetxt(f"{out}/nwater.dat", np.array(nwat), header="time(ns) n_water_within_0.35nm", fmt="%10.3f")
ins = dg[dg < R90]
# A solute that never reached the core leaves the guest histogram empty. Reading
# a peak position off it returns the first bin, and the outer-half fraction
# divides by nothing, so the report read as a broken calculation rather than as
# a solute that simply stayed outside. Those two lines are only written when
# there is something inside to describe.
inside = ins.size
with open(f"{out}/localisation.txt", "w") as f:
    f.write(f"frames                 {n}\nsamples                {dg.size}\n")
    f.write(f"R_core (90%ile PPO)    {R90:.3f} nm\n")
    f.write(f"<d>/R_core             {dg.mean()/R90:.3f}   (uniform sphere = 0.750)\n")
    f.write(f"fraction inside R_core {100 * inside / dg.size:.1f} %\n")
    if inside:
        f.write(f"peak of number density {cent[np.argmax(acc['guest'])]:.3f} nm\n")
        f.write(f"  of those, in outer half-volume "
                f"{100*(ins>R90/2**(1/3)).mean():.1f} %  (uniform = 50 %)\n")
    else:
        f.write(f"the solute never came within R_core over these {n} frames. "
                f"Its closest approach was {dg.min():.2f} nm against an "
                f"R_core of {R90:.2f} nm.\n")
        f.write("This is the expected reading for a short run of the solution "
                "route; uptake in the paper took tens of nanoseconds.\n")
print(open(f"{out}/localisation.txt").read())
