Nocc    = 137
occAct  = [131,135,136] # Put the numbers here!
virAct  = [138,137] # Put the numbers here!
allOccs = list(range(Nocc))
frozen  = [] # just leave it empty
doccs   = list(filter(lambda x: (occAct.count(x)==0), allOccs))
acts    = occAct + virAct
