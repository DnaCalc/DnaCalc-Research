from xl import *
rows=[]
b=bits(1.3)
for i in range(20000):
    rows.append([unb(b+i)])
run("GAMMA",rows,"q21.jsonl",chunk=6000)
