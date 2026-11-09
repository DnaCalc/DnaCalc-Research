def norm(e):
    if e is None: return e
    if "#NUM!" in e: return "error:Num"
    if "#N/A" in e: return "error:NA"
    if "#DIV/0!" in e: return "error:Div0"
    if "#VALUE!" in e: return "error:Value"
    return e
