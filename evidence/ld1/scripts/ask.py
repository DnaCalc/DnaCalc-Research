import sys, json, struct
sys.path.insert(0, "/home/"+__import__("getpass").getuser()+"/projects/DnaCalc/DnaCalc/tools/excel-service")
from client import Client, number, text, boolean, blank, error, number_array
def h(x): return struct.pack(">d", x).hex()
def f(hx): return struct.unpack(">d", bytes.fromhex(hx))[0]
excel = Client()
def ask(fn, rows, **kw):
    kw.setdefault("priority", "interactive" if len(rows) <= 100 else "bulk")
    r = excel.ask(fn, [[number(v) if isinstance(v,(int,float)) else v for v in row] for row in rows], **kw)
    return r
