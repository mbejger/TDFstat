import numpy as np
import matplotlib.pyplot as plt 
#import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter,MultipleLocator,AutoMinorLocator
from matplotlib.pyplot import figure
import matplotlib.patches as patches

import os 
import sys 

# +- freq. bins - for plotting 
fshift = int(sys.argv[2])  

chunk_size = 1000000

 
def read_triggers(filename, fshift): 

    rectype = np.dtype( [ ('f', 'float32'), ('fdot', 'float32'), 
                          ('fa', 'float32'), ('fda', 'float32'),
                          ('snr', 'float32') ] )

    with open(filename) as file:

        print("Filename: {}".format(filename)) 
        offset = 0 
        while True:
            recs = None
            file.seek(offset*rectype.itemsize)
            recs = np.fromfile(file, dtype=rectype, count=chunk_size)
            offset += chunk_size
            
            if len(recs) < chunk_size:
                break

    fstat = np.split(np.asarray([r[4] for r in recs]), int(2*fshift+1))
    max_fstat = np.max(fstat)   

    out = 1 - fstat/max_fstat
 
    return out, max_fstat 



if os.stat(sys.argv[1]).st_size == 0:
    print("{} size is zero. Exiting...".format(sys.argv[1])) 
    sys.exit() 


out, max_fstat = read_triggers(sys.argv[1], fshift)
out_rgb = np.stack([out]*3, axis=-1)
out_rgb[fshift, len(out)] = [1, 0, 0]  # mark the center in red

fig = plt.figure(constrained_layout=True)
ax = fig.subplot_mosaic("""a""")

ax["a"].set_adjustable("box")
p = ax["a"].imshow(out_rgb, interpolation='none')

ax["a"].set_title(sys.argv[1] + " | max Fstat: " + str(max_fstat), fontsize=8) 

if ("-show" in sys.argv[2:] ):
    plt.show()
else:
    try:
        #plt.savefig(sys.argv[-1] +".png", format="png", dpi=150, bbox_inches='tight')
        plt.savefig(sys.argv[-1] +".pdf", bbox_inches='tight')
    except:
        print("Warning: can't save ", sys.argv[-1] + ".png")

