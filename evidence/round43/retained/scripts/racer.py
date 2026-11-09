import os, subprocess, tempfile, struct
RACER=os.path.expanduser("~/codex-runs/r43/racer-base")
def f(b): return struct.unpack("<d", struct.pack("<Q", int(b) & 0xffffffffffffffff))[0]
def h16(x): return "%016x"%struct.unpack("<Q",struct.pack("<d",float(x)))[0]
def invoke(rows):
    d=tempfile.mkdtemp(dir=os.path.expanduser("~/codex-runs/r43/tmp"))
    fi=os.path.join(d,"in"); fo=os.path.join(d,"out")
    with open(fi,"w") as fh:
        for fn,args in rows: fh.write(fn+" "+" ".join(a if isinstance(a,str) else h16(a) for a in args)+"\n")
    subprocess.run([RACER,"invoke-rows",fi,fo],check=True,capture_output=True)
    out=[]
    for line in open(fo):
        s=line.strip(); out.append(f(int(s,16)) if s.startswith("0x") else s)
    os.remove(fi); os.remove(fo); os.rmdir(d); return out
