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
#                print('file offset = ',offset*rectype.itemsize,' npoints = ',len(recs['f']), end='\r')
#                print('\nend of file')
                break

    snr = np.split(np.asarray([r[4] for r in recs]), int(2*fshift+1))
    
    return snr 



if os.stat(sys.argv[1]).st_size == 0:
    print("{} size is zero. Exiting...".format(sys.argv[1])) 
    sys.exit() 


snr = read_triggers(sys.argv[1], fshift)


fig = plt.figure(constrained_layout=True)
ax = fig.subplot_mosaic("""a""")

ax["a"].set_adjustable("box")

p = ax["a"].imshow(snr, cmap='binary')

## adding injection point 
#inj = ax["a"].scatter(len(snr)-1, fshift, s=5.5, c='red', marker='s')

#fig.colorbar(p, ax=ax["a"])

#ax["a"].set_xticks([])
#ax["a"].set_yticks([])

#ax["a"].axis("off") 

# Create a Rectangle patch
#rect = patches.Rectangle((fshift-fcellsize/2, fshift-fdotcellsize/2), fcellsize, fdotcellsize, linewidth=1, edgecolor='g', facecolor='none')

# Add the patch to the Axes
#ax["a"].add_patch(rect)

ax["a"].set_title(sys.argv[1]) 

if ("-show" in sys.argv[2:] ):
    plt.show()
else:
    try:
        plt.savefig(sys.argv[-1] +".png", format="png", dpi=150, bbox_inches='tight')
    except:
        print("Warning: can't save ", sys.argv[-1] + ".png")

