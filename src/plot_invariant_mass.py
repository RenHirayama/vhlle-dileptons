import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input", "-i", nargs='+',required=True)
args = parser.parse_args()

def midpoints(x):
    return (x[1:]+x[:-1])/2

mass_bins=np.linspace(0,2.5,100)
for f in args.input:
    fulldata = np.loadtxt(f)#,usecols=[0,1,2,3],dtype=[('ityp',np.int16),('weight',np.float64),
                        #                        ('mass',np.float64),('p0_el',np.float64)])

    hist, bins = np.histogram(fulldata[:,2],weights=fulldata[:,1],bins=mass_bins)
    plt.plot(midpoints(bins),hist,'-k',label='sum',lw=2.5)
#    data = data[np.where(data[:,0]==103)]
    ptypes = np.unique(fulldata[:,0])
    for ityp in ptypes:
        data = fulldata[np.where(fulldata[:,0]==ityp)]
        hist, bins = np.histogram(data[:,2],weights=data[:,1],bins=mass_bins)
        plt.plot(midpoints(bins),hist,label=int(ityp))
plt.legend(loc=0)
plt.yscale('log')
plt.savefig('nodelta.png',dpi=500)
plt.show()