from xl import *
import sys
for x0,tag in ((1.301,"a"),(1.302,"b"),(1.304,"c"),(1.308,"d"),(1.316,"e"),(1.332,"f"),(1.364,"g")):
    rows=[]
    b=bits(x0)
    for i in range(6000): rows.append([unb(b+i)])
    run("GAMMA",rows,"q24%s.jsonl"%tag,chunk=6000)
