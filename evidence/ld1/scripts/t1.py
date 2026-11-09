from ask import *
r = ask("CUMPRINC", [[0.05/12, 360, 200000, 1, 12, 0],[0.05/12,360,200000,5,5,0]])
print(r["excel_build"], r.get("env_id"), r.get("excel_version"))
for row in r["rows"]: print(row)
r = ask("PPMT", [[0.05/12, 5, 360, 200000, 0, 0]])
print(r["rows"])
