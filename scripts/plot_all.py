"""Draw the .dat and .xvg files in a results directory.  plot_all.py OUTDIR"""
import sys, os, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
O = sys.argv[1]; F = os.path.join(O,'figures'); os.makedirs(F, exist_ok=True)
NAVY='#2c3e50'; RED='#c0392b'; BLUE='#2980b9'; GREY='#95a5a6'
plt.rcParams.update({'font.size':13,'axes.linewidth':1.6,'font.weight':'bold',
                     'axes.labelweight':'bold','axes.labelsize':15})
def xvg(f):
    t,v=[],[]
    for L in open(f,errors='ignore'):
        if L[:1] in '#@': continue
        p=L.split()
        try: t.append(float(p[0])); v.append([float(x) for x in p[1:]])
        except: pass
    return np.array(t), np.array(v)
def table(f):
    """loadtxt returns a 1-D array for a single data row, so it is made 2-D."""
    d = np.loadtxt(f)
    d = np.atleast_2d(d)
    return d if d.size else None
def fin(ax):
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5)); ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    for t in ax.get_xticklabels()+ax.get_yticklabels(): t.set_fontweight('bold')
def key(ax, n):
    """A legend below the axes, where it cannot sit on a curve."""
    ax.legend(frameon=False, ncol=n, loc='upper center', bbox_to_anchor=(0.5, -0.18),
              fontsize=11, handlelength=1.8, columnspacing=1.4)

p=os.path.join(O,'radial_density.dat')
if os.path.exists(p) and table(p) is not None:
    d=table(p); fig,ax=plt.subplots(figsize=(7,5))
    for i,(lab,c) in enumerate([('core',NAVY),('corona',BLUE),('solute',RED),('water',GREY)],start=1):
        y=d[:,i]; ax.plot(d[:,0], y/max(y.max(),1e-9), color=c, lw=3, label=lab)
    ax.set_xlabel('r from micelle COM (nm)'); ax.set_ylabel('normalised density')
    fin(ax); key(ax, 4); fig.tight_layout(); fig.savefig(f"{F}/radial_density.png",dpi=200)
for nm,xl,yl in [('encapsulation','time (ns)','solutes inside R$_{core}$'),
                 ('nwater','time (ns)','water within 0.35 nm')]:
    p=os.path.join(O,f"{nm}.dat")
    if os.path.exists(p):
        d=table(p)
        if d is None: continue
        fig,ax=plt.subplots(figsize=(7,4.5))
        ax.plot(d[:,0],d[:,1],color=RED,lw=2.2,marker='o' if len(d)<3 else None); ax.set_xlabel(xl); ax.set_ylabel(yl)
        fin(ax)
        if nm=='encapsulation':
            # a count, so the axis carries whole numbers. Must come after fin(),
            # which sets its own locator.
            top=max(int(d[:,1].max()), 1)
            ax.set_ylim(-0.04*top, top*1.12)
            ax.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=min(top+1,6)))
        fig.tight_layout(); fig.savefig(f"{F}/{nm}.png",dpi=200)
for nm,xl in [('assoc_distance','centroid distance (nm)'),('assoc_angle','interplanar angle (deg)')]:
    p=os.path.join(O,f"{nm}.dat")
    if os.path.exists(p):
        d=table(p)
        if d is None: continue
        fig,ax=plt.subplots(figsize=(7,4.5))
        ax.plot(d[:,0],d[:,2],color=BLUE,lw=3); ax.fill_between(d[:,0],0,d[:,2],color=BLUE,alpha=.25)
        ax.set_xlabel(xl); ax.set_ylabel('normalised density')
        fin(ax); fig.tight_layout(); fig.savefig(f"{F}/{nm}.png",dpi=200)
for nm,yl in [('rg','R$_g$ (nm)'),('sasa','SASA (nm$^2$)')]:
    p=os.path.join(O,f"{nm}.xvg")
    if os.path.exists(p):
        t,v=xvg(p)
        if not v.size: continue
        fig,ax=plt.subplots(figsize=(7,4.5))
        ax.plot(t/1000,v[:,0],marker='o' if len(t)<3 else None,color=NAVY,lw=2); ax.set_xlabel('time (ns)'); ax.set_ylabel(yl)
        fin(ax); fig.tight_layout(); fig.savefig(f"{F}/{nm}.png",dpi=200)
p=os.path.join(O,'rdf.xvg')
if os.path.exists(p) and xvg(p)[1].size:
    t,v=xvg(p)
    fig,ax=plt.subplots(figsize=(7,5))
    for i,(lab,c) in enumerate([('core',NAVY),('corona',BLUE),('water',GREY)]):
        if i<v.shape[1]: ax.plot(t,v[:,i],color=c,lw=3,label=lab)
    ax.set_xlabel('r (nm)'); ax.set_ylabel('g(r)')
    fin(ax); key(ax, 3); fig.tight_layout(); fig.savefig(f"{F}/rdf.png",dpi=200)
p=os.path.join(O,'lie.xvg')
if os.path.exists(p) and xvg(p)[1].size:
    t,v=xvg(p)
    fig,ax=plt.subplots(figsize=(7.5,5))
    lab=['Coul: solute-core','vdW: solute-core','Coul: solute-water','vdW: solute-water']
    for i in range(min(4,v.shape[1])):
        ax.plot(t/1000,v[:,i],lw=2,label=lab[i],marker='o' if len(t)<3 else None,
                color=[RED,NAVY,'#e67e22',GREY][i], ls='--' if i%2==0 else '-')
    ax.set_xlabel('time (ns)'); ax.set_ylabel('interaction energy (kJ mol$^{-1}$)')
    fin(ax); key(ax, 2); fig.tight_layout(); fig.savefig(f"{F}/lie.png",dpi=200)
print("  figures:", sorted(os.listdir(F)))
