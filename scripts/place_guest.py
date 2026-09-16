"""Place solute molecules inside the micelle (shell) or in bulk water (solution).
usage: place_guest.py host.gro solute.pdb N {shell|soak} out.gro [seed] [box_nm]

box_nm is the edge of the cubic box the system will end up in. Without it the
placement radius comes from the host's own box."""
import sys, numpy as np

def read_gro(fn):
    L=open(fn).read().splitlines(); n=int(L[1])
    rec=[(l[:5],l[5:10],l[10:15],l[15:20],
          float(l[20:28]),float(l[28:36]),float(l[36:44])) for l in L[2:2+n]]
    box=np.array([float(x) for x in L[2+n].split()[:3]])
    return L[0], rec, box

def read_guest(fn):
    xs, names = [], []
    for l in open(fn):
        if l.startswith(('ATOM','HETATM')):
            xs.append((float(l[30:38]),float(l[38:46]),float(l[46:54])))
            names.append((l[12:16].strip(), l[17:20].strip()))
        elif fn.endswith('.gro'):
            pass
    if not xs:                                   # a gro was given instead of a pdb
        L=open(fn).read().splitlines(); n=int(L[1])
        for l in L[2:2+n]:
            xs.append((float(l[20:28])*10,float(l[28:36])*10,float(l[36:44])*10))
            names.append((l[10:15].strip(), l[5:10].strip()))
    return np.array(xs)/10.0, names               # nm

def rand_rot(rng):
    q=rng.normal(size=4); q/=np.linalg.norm(q); w,x,y,z=q
    return np.array([[1-2*(y*y+z*z),2*(x*y-z*w),2*(x*z+y*w)],
                     [2*(x*y+z*w),1-2*(x*x+z*z),2*(y*z-x*w)],
                     [2*(x*z-y*w),2*(y*z+x*w),1-2*(x*x+y*y)]])

def main():
    host_f, guest_f, N, mode, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5]
    seed = int(sys.argv[6]) if len(sys.argv)>6 else 2026
    rng = np.random.default_rng(seed)
    title, rec, box = read_gro(host_f)
    host_box = box.copy()                      # the box the host was equilibrated in
    if len(sys.argv) > 7:                      # the box the system will end up in
        box = np.full(3, float(sys.argv[7]))
    hp = np.array([[r[4],r[5],r[6]] for r in rec])
    g0, gnames = read_guest(guest_f)
    g0 = g0 - g0.mean(0)
    resn = gnames[0][1] if gnames[0][1] else 'LIG'

    com = hp.mean(0)
    d = np.linalg.norm(hp-com, axis=1)
    R95 = np.percentile(d, 95)

    # The test is against the box the host was equilibrated in, not the furthest
    # atom: an equilibrated corona reaching past half its own box is normal.
    if box.min() < host_box.min() - 0.05:
        out_n = int((d > box.min() / 2).sum())
        print(f"  [warn] the requested {box.min():.1f} nm box is smaller than the"
              f" {host_box.min():.1f} nm the host was equilibrated in, so"
              f" {out_n} atoms reach past half the box and wrap onto the far"
              f" side.")
        print(f"         The corona then interdigitates with its own periodic"
              f" image, which is how a concentrated system is meant to look but"
              f" is not what a dilute one should do. Check that the box after"
              f" the pressure step is close to {box.min():.1f} nm. To start from"
              f" a host that is already this dense, compress one with"
              f" densify.py.")
    if mode == 'shell':
        # the template's hollow is 4.6 nm across, so a solute centre stays
        # inside 1.5 nm and leaves room for the molecule itself
        rmax = 1.5; lo, hi = 0.0, rmax
        print(f"  shell: centres within {hi:.1f} nm of the middle")
    elif mode == 'soak':
        # anywhere in the box that is far enough from the polymer, corners
        # included
        lo = hi = None
        print(f"  solution: anywhere in the water (micelle R95 {R95:.1f} nm, box {box.min():.1f} nm)")
    else: sys.exit("mode must be shell or soak")

    # The host never moves, so the clash test is prepared once, keeping only
    # host atoms that can reach the placement region.
    def _min_dist(ref, xyz, box):
        dd = ref - xyz[:, None, :]
        dd -= box * np.round(dd / box)
        return np.linalg.norm(dd, axis=2).min()

    try:
        from scipy.spatial import cKDTree
    except ImportError:
        cKDTree = None
    reach = np.linalg.norm(g0, axis=1).max() + 0.5
    if lo is None:
        host_near = hp                                   # soak: anywhere in the box
    else:
        dh = hp - com
        dh -= box * np.round(dh / box)
        host_near = hp[np.linalg.norm(dh, axis=1) <= hi + reach + 0.5]
    # boxsize makes the tree periodic, so a molecule near a face still sees its
    # image. It needs coordinates inside [0, L), which is where a gro keeps them.
    tree = None
    if cKDTree is not None:
        try:
            tree = cKDTree(np.mod(host_near, box), boxsize=box)
        except Exception:
            tree = None
    placed, placed_pts = [], np.zeros((0, 3))
    tries = 0
    while len(placed) < N and tries < 200000:
        tries += 1
        if lo is None:                      # soak, uniform over the box
            c = com + (rng.random(3) - 0.5) * (box - 1.0)
        else:                               # shell, uniform in a sphere
            u = rng.normal(size=3); u /= np.linalg.norm(u)
            c = com + u * (lo**3 + rng.random()*(hi**3-lo**3))**(1/3)
        xyz = g0 @ rand_rot(rng).T + c
        # the neighbour tree answers the clash test without a full distance array
        if tree is not None:
            if tree.query_ball_point(xyz, 0.28, return_length=True).any():
                continue
            if placed_pts.size and _min_dist(placed_pts, xyz, box) < 0.28:
                continue
            placed.append(xyz)
            placed_pts = np.vstack((placed_pts, xyz))
            continue
        ref = host_near if placed_pts.size == 0 else np.vstack((host_near, placed_pts))
        dd = ref - xyz[:,None,:]
        dd -= box*np.round(dd/box)
        # 2.8 A against anything already there. For the solution route also keep
        # clear of the polymer surface, so the solute starts in water and is not
        # already touching the micelle it is supposed to find on its own.
        near = np.linalg.norm(dd, axis=2).min()
        if near < 0.28:
            continue
        if lo is None:
            dh = hp - xyz[:, None, :]
            dh -= box * np.round(dh / box)
            if np.linalg.norm(dh, axis=2).min() < 0.60:
                continue
        placed.append(xyz)
        placed_pts = np.vstack((placed_pts, xyz))
    if len(placed) < N:
        sys.exit(f"  [failed] placed only {len(placed)} of {N}. Use fewer, or a larger box.")

    lines, rid, aid = [], 0, 0
    for r in rec:
        if r[0].strip() and (not lines or r[0]!=rec[0][0]): pass
    prev=None
    for r in rec:
        if r[0] != prev: rid += 1; prev = r[0]
        aid += 1
        lines.append(f"{rid%100000:5d}{r[1]}{r[2]}{aid%100000:5d}{r[4]:8.3f}{r[5]:8.3f}{r[6]:8.3f}")
    for m in placed:
        rid += 1
        for (an,_),p in zip(gnames, m):
            aid += 1
            lines.append(f"{rid%100000:5d}{resn:>5s}{an:>5s}{aid%100000:5d}{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}")
    with open(out,'w') as f:
        f.write(f"{title} + {N} {resn}\n{len(lines)}\n")
        f.write("\n".join(lines)+"\n")
        f.write(f"{box[0]:10.5f}{box[1]:10.5f}{box[2]:10.5f}\n")
    print(f"  -> {out}  ({len(lines)} atoms, {N} x {resn}, {tries} tries)")
    print(f"RESNAME={resn}")

if __name__ == '__main__': main()
