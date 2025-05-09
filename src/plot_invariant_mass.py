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

mass_bins=np.linspace(0,2.5,100)
x_bins = midpoints(mass_bins)
for i_file,f in enumerate(args.input):
    fulldata = np.loadtxt(f)
    hist, bins = np.histogram(fulldata[:,3],weights=fulldata[:,2],bins=mass_bins)
    plt.plot(x_bins, hist/(2*x_bins), '--',label='sum '+f.split('_')[1].split('f')[0],lw=2.5)
    print(fulldata)
    ptypes = np.unique(fulldata[:,1]).astype(int)
    print(ptypes)
    for ityp in ptypes:
        data = fulldata[np.where(fulldata[:,1].astype(int)==ityp)]
        print(data[data[:,2]==0])
        hist, bins = np.histogram(data[:,3],weights=data[:,2],bins=mass_bins)
        plt.plot(x_bins,hist/(2*x_bins),label=f.split('_')[1].split('f')[0]+" "+str(int(ityp)))
plt.legend(loc=0)
plt.xlim(0,1.4)
plt.yscale('log')
#plt.ylim(10**-12,10**-1)
plt.savefig('diltest_3.5.png',dpi=500)
plt.show()