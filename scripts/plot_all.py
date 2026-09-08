"""결과 폴더의 .dat/.xvg 를 읽어 PNG 를 만든다.  plot_all.py OUTDIR"""
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
def fin(ax):
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5)); ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    for t in ax.get_xticklabels()+ax.get_yticklabels(): t.set_fontweight('bold')

p=os.path.join(O,'radial_density.dat')
if os.path.exists(p):
    d=np.loadtxt(p); fig,ax=plt.subplots(figsize=(7,5))
    for i,(lab,c) in enumerate([('core',NAVY),('corona',BLUE),('guest',RED),('water',GREY)],start=1):
        y=d[:,i]; ax.plot(d[:,0], y/max(y.max(),1e-9), color=c, lw=3, label=lab)
    ax.set_xlabel('r from micelle COM (nm)'); ax.set_ylabel('normalised density')
    ax.legend(frameon=False); fin(ax); fig.tight_layout(); fig.savefig(f"{F}/radial_density.png",dpi=200)
for nm,xl,yl in [('encapsulation','time (ns)','guests inside R$_{core}$'),
                 ('nwater','time (ns)','water within 0.35 nm')]:
    p=os.path.join(O,f"{nm}.dat")
    if os.path.exists(p):
        d=np.loadtxt(p); fig,ax=plt.subplots(figsize=(7,4.5))
        ax.plot(d[:,0],d[:,1],color=RED,lw=2.2); ax.set_xlabel(xl); ax.set_ylabel(yl)
        fin(ax); fig.tight_layout(); fig.savefig(f"{F}/{nm}.png",dpi=200)
for nm,xl in [('assoc_distance','centroid distance (nm)'),('assoc_angle','interplanar angle (deg)')]:
    p=os.path.join(O,f"{nm}.dat")
    if os.path.exists(p):
        d=np.loadtxt(p); fig,ax=plt.subplots(figsize=(7,4.5))
        ax.plot(d[:,0],d[:,2],color=BLUE,lw=3); ax.fill_between(d[:,0],0,d[:,2],color=BLUE,alpha=.25)
        ax.set_xlabel(xl); ax.set_ylabel('normalised density')
        fin(ax); fig.tight_layout(); fig.savefig(f"{F}/{nm}.png",dpi=200)
for nm,yl in [('rg','R$_g$ (nm)'),('sasa','SASA (nm$^2$)')]:
    p=os.path.join(O,f"{nm}.xvg")
    if os.path.exists(p):
        t,v=xvg(p); fig,ax=plt.subplots(figsize=(7,4.5))
        ax.plot(t/1000,v[:,0],color=NAVY,lw=2); ax.set_xlabel('time (ns)'); ax.set_ylabel(yl)
        fin(ax); fig.tight_layout(); fig.savefig(f"{F}/{nm}.png",dpi=200)
p=os.path.join(O,'rdf.xvg')
if os.path.exists(p):
    t,v=xvg(p); fig,ax=plt.subplots(figsize=(7,5))
    for i,(lab,c) in enumerate([('core',NAVY),('corona',BLUE),('water',GREY)]):
        if i<v.shape[1]: ax.plot(t,v[:,i],color=c,lw=3,label=lab)
    ax.set_xlabel('r (nm)'); ax.set_ylabel('g(r)'); ax.legend(frameon=False)
    fin(ax); fig.tight_layout(); fig.savefig(f"{F}/rdf.png",dpi=200)
p=os.path.join(O,'lie.xvg')
if os.path.exists(p):
    t,v=xvg(p); fig,ax=plt.subplots(figsize=(7.5,5))
    lab=['Coul: guest-core','vdW: guest-core','Coul: guest-water','vdW: guest-water']
    for i in range(min(4,v.shape[1])):
        ax.plot(t/1000,v[:,i],lw=2,label=lab[i],
                color=[RED,NAVY,'#e67e22',GREY][i], ls='--' if i%2==0 else '-')
    ax.set_xlabel('time (ns)'); ax.set_ylabel('interaction energy (kJ mol$^{-1}$)')
    ax.legend(frameon=False,fontsize=11); fin(ax); fig.tight_layout(); fig.savefig(f"{F}/lie.png",dpi=200)
print("  그림:", sorted(os.listdir(F)))
