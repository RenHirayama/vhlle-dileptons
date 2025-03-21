import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input", "-i", nargs='+',required=True)
args = parser.parse_args()

def midpoints(x):
    return (x[1:]+x[:-1])/2

def pT(px,py):
    return np.sqrt(px**2+py**2)

mass_bins=np.linspace(0,5,100)
x_bins = midpoints(mass_bins)
for i_file,f in enumerate(args.input):
    fulldata = np.loadtxt(f)
    fulldata = fulldata[np.where(pT(fulldata[:,4],fulldata[:,5])>0)[0]]
    hist, bins = np.histogram(pT(fulldata[:,4],fulldata[:,5]),weights=fulldata[:,1],bins=mass_bins)
    plt.plot(x_bins, hist/(2*x_bins), '--',label='sum '+f.split('_')[1].split('f')[0],lw=2.5)

    ptypes = np.unique(fulldata[:,0])
    for ityp in ptypes:
        data = fulldata[np.where(fulldata[:,0]==ityp)]
        hist, bins = np.histogram(pT(data[:,4],data[:,5]),weights=data[:,1],bins=mass_bins)
        plt.plot(x_bins,hist/(2*x_bins),label=f.split('_')[1].split('f')[0]+" "+str(int(ityp)))
plt.legend(loc=0)
plt.xlim(0,2)
plt.yscale('log')
plt.ylim(10**-12,10**-1)
plt.savefig('mupion0_coldest2.png',dpi=500)
plt.show()