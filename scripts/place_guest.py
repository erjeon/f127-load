"""화물 분자를 미셀 안(shell) 또는 벌크 물 영역(soak)에 무작위 배치한다.
사용: place_guest.py host.gro guest.pdb N {shell|soak} out.gro [seed]
겹침 검사를 통과할 때까지 위치·회전을 다시 뽑는다."""
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
    if not xs:                                   # gro 입력 지원
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
    hp = np.array([[r[4],r[5],r[6]] for r in rec])
    g0, gnames = read_guest(guest_f)
    g0 = g0 - g0.mean(0)
    resn = gnames[0][1] if gnames[0][1] else 'LIG'

    com = hp.mean(0)
    d = np.linalg.norm(hp-com, axis=1)
    R95 = np.percentile(d, 95)
    if mode == 'shell':
        rmax = 1.5; lo, hi = 0.0, rmax
        print(f"  shell 배치: COM 반경 {hi:.1f} nm 이내")
    elif mode == 'soak':
        lo, hi = R95 + 1.0, box.min()/2 - 0.8
        if lo >= hi: sys.exit(f"  [FAIL] 벌크 영역이 없다 (미셀 R95={R95:.1f}, 박스 {box.min():.1f} nm). 박스를 키워라.")
        print(f"  soak 배치: COM 반경 {lo:.1f}–{hi:.1f} nm")
    else: sys.exit("mode 는 shell 또는 soak")

    placed, occ = [], list(hp)
    tries = 0
    while len(placed) < N and tries < 200000:
        tries += 1
        u = rng.normal(size=3); u /= np.linalg.norm(u)
        r = (lo**3 + rng.random()*(hi**3-lo**3))**(1/3)
        c = com + u*r
        xyz = g0 @ rand_rot(rng).T + c
        ref = np.array(occ)
        dd = ref - xyz[:,None,:]
        dd -= box*np.round(dd/box)
        if np.linalg.norm(dd,axis=2).min() < 0.28:   # 2.8 A 이내면 재시도
            continue
        placed.append(xyz); occ.extend(list(xyz))
    if len(placed) < N:
        sys.exit(f"  [FAIL] {N}개 중 {len(placed)}개만 배치. 반경/개수를 줄여라.")

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
    print(f"  → {out}  ({len(lines)} atoms, {N} x {resn}, {tries} tries)")
    print(f"RESNAME={resn}")

if __name__ == '__main__': main()
