"""물/이온을 무작위로 제거해 목표 wt% 로 만든다. 박스는 그대로 두고 NPT 가 압축하게 한다."""
import sys, random
src, dst, top, wt = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
random.seed(20260903)
NMIC = 70010                      # S1P1x34 + PYRx8
MW_POL = 34*12586.0; MW_W = 18.015; NA = 6.02214076e23

L = open(src).read().splitlines()
n = int(L[1]); head, atoms, box = L[0], L[2:2+n], L[2+n]

mic = atoms[:NMIC]
rest = atoms[NMIC:]
# 잔기 단위로 묶기
res, cur, curid = [], [], None
for a in rest:
    rid = a[:5]
    if rid != curid and cur: res.append(cur); cur=[]
    curid = rid; cur.append(a)
if cur: res.append(cur)
wat = [r for r in res if r[0][5:10].strip()=='TIP3']
sod = [r for r in res if r[0][5:10].strip()=='SOD']
cla = [r for r in res if r[0][5:10].strip()=='CLA']
print(f"  원본: TIP3 {len(wat)}  SOD {len(sod)}  CLA {len(cla)}")

m_pol = MW_POL/NA
n_wat = int(round(m_pol*(1-wt)/wt / (MW_W/NA)))
V_L   = (m_pol/wt)/1.0 * 1e-3            # g / (g cm-3) = cm3 -> L
n_ion = int(round(0.154 * V_L * NA))
if n_wat > len(wat): sys.exit(f"  [FAIL] 물이 부족하다: 필요 {n_wat} > 보유 {len(wat)}")
print(f"  목표 wt%={wt*100:.0f} → TIP3 {n_wat}  SOD/CLA {n_ion}  (예상 부피 {V_L*1e21:.0f} nm3, box {(V_L*1e21)**(1/3):.2f} nm)")

keep_w = sorted(random.sample(range(len(wat)), n_wat))
keep_i = sorted(random.sample(range(len(sod)), min(n_ion, len(sod))))
keep_c = sorted(random.sample(range(len(cla)), min(n_ion, len(cla))))
out = list(mic)
for i in keep_w: out += wat[i]
for i in keep_i: out += sod[i]
for i in keep_c: out += cla[i]

# 잔기/원자 번호 재부여
lines, rid, aid = [], 0, 0
prev = None
for a in out:
    if a[:5] != prev: rid += 1; prev = a[:5]
    aid += 1
    lines.append(f"{rid%100000:5d}{a[5:15]}{aid%100000:5d}{a[20:]}")
with open(dst,'w') as f:
    f.write(head+"\n"); f.write(f"{len(lines)}\n")
    f.write("\n".join(lines)+"\n"); f.write(box+"\n")

t = open(top).read().splitlines()
o = []
for L2 in t:
    k = L2.split()
    if len(k)==2 and k[0] in ('TIP3','SOD','CLA'):
        o.append(f"{k[0]:<12s}{ {'TIP3':n_wat,'SOD':len(keep_i),'CLA':len(keep_c)}[k[0]] }")
    else: o.append(L2)
open(top,'w').write("\n".join(o)+"\n")
print(f"  → {dst}  총 {len(lines)} atoms")
