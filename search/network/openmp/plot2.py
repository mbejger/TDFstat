import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.pyplot import figure

import os 
import sys 


def read_triggers(filename): 

    # f, s(pindown), d(eclination), r(ight ascension), v(alue of Fstat) 
    rectype = np.dtype( [ ('f', 'float32'), ('s', 'float32'), 
                          ('d', 'float32'), ('r', 'float32'),
                          ('v', 'float32') ] )


    chunk_size = 1048576
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


    sky_pos = np.unique(recs[['d', 'r']])

    out = []
    fstat_max = 0 
    for pos in sky_pos:

            d, r = pos

            # Filter (f, fdot, fstat) with the same (de, ra) values
            fsv = recs[(recs['d'] == d) & (recs['r'] == r)][['f', 's', 'v']]
            fdot = np.unique(fsv[['s']])

            fstat = [] 
            for s in fdot:
                v = np.array(recs[(recs['d'] == d) & (recs['r'] == r) & (recs['s'] == s[0])][['v']], np.dtype('float32'))
 
                fstat.append(v) 

                vmax = max(v) 
            
                if vmax > fstat_max: 
                    fstat_max = vmax 

            out.append([pos, fstat])
    
    return out, fstat_max 


if os.stat(sys.argv[1]).st_size == 0:
    print("{} size is zero. Exiting...".format(sys.argv[1])) 
    sys.exit() 


out, fstat_max = read_triggers(sys.argv[1])


print(fstat_max) 

#sys.exit()

gsize = 2
skypos = 2*gsize + 1

# currently first 9 sky positions 
#fig, ax = plt.subplots(skypos, skypos)
#plt.subplots_adjust(top = 0.99, bottom=0.01, hspace=0, wspace=-15) 
#fig.tight_layout(pad=-10)

fig = plt.figure()
gs = fig.add_gridspec(skypos, skypos, hspace=0, wspace=-0.1)
ax = gs.subplots(sharex='col', sharey='row')

for i in range(skypos):
    for j in range(skypos):

        snr = out[i*skypos+j][1]/fstat_max 
        maxsnr = np.max(snr) 
        print(out[i*skypos+j][0]) 
 
        ax[i,j].imshow(snr, cmap='binary', vmin=0, vmax=1)
        ax[i,j].set_xticks([])
        ax[i,j].set_yticks([])

        ax[i,j].text(10, 10, "{:.2f}".format(float(maxsnr)))

#if ("-show" in sys.argv[2:] ):
#    plt.show()
#else:
#    try:

plt.savefig(sys.argv[-1] +".pdf", format="pdf", bbox_inches='tight')
#    except:
#        print("Warning: can't save ", sys.argv[-1] + ".png")

