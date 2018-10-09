import sys
import numpy as np

fn = sys.argv[1]
if fn == "stdin":
    fhandle = sys.stdin
else:
    fhandle = open(sys.argv[1])

for f in fhandle:
    if "n" in f:
        for line in fhandle:
            N = int(line.split()[0])
            sfs = np.zeros([N,N,N])
            break
    elif "//" in f:
        print(sfs.flatten())
    else:
        toks = f.strip().split()
        sfs[int(toks[1]),int(toks[2]),int(toks[3])] += 1



