import sys, MDAnalysis as mda
gro, out = sys.argv[1], sys.argv[2]
u = mda.Universe(gro)
g = {
 'micelle': u.select_atoms("resname PROXS PROXR PROX ETHOX ETHO"),
 'PYR'    : u.select_atoms("resname PYR"),
 'W_ION'  : u.select_atoms("resname TIP3 SOL SOD CLA POT"),
}
tot = sum(len(v) for v in g.values())
assert tot == len(u.atoms), f"group sum {tot} != {len(u.atoms)}"
with open(out,'w') as f:
    f.write("[ System ]\n")
    for i,a in enumerate(u.atoms.indices+1):
        f.write(f"{a:6d}" + ("\n" if (i+1)%15==0 else " "))
    f.write("\n")
    for k,v in g.items():
        f.write(f"[ {k} ]\n")
        for i,a in enumerate(v.indices+1):
            f.write(f"{a:6d}" + ("\n" if (i+1)%15==0 else " "))
        f.write("\n")
print("ndx ok:", {k:len(v) for k,v in g.items()})
