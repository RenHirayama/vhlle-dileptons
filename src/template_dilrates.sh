#!/bin/bash
#SBATCH --job-name=dil_vhlle_ENERGY-EVNUM
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --output=gen_dil_vhlle_ENERGY-EVNUM-%j.out
#SBATCH --error=gen_dil_vhlle_ENERGY-EVNUM-%j.err
#SBATCH --partition=long
#SBATCH --time=3-5:29:00
#SBATCH --mail-user=hirayama@itp.uni-frankfurt.de
#SBATCH --mail-type=FAIL
#SBATCH --mem=12GB

SOURCE="/lustre/hyihp/rhirayam/smash-vhlle-hybrid/vhlle-dileptons/src"

export model=2
export ecm=ENERGY
export betalab=0
export vacuum=0
export eos=6
export rates=4
export pikchem=1
export baryons=1
# WRITE DILEPTON / PHOTON OUTPUT f71? *1* = YES, *0* = NO
# (IF NO, ONLY CELL PROPERTIES AND NON-THERMAL f72 CALCULATED)
export output=1
export leptype=0
export fourpimix=1
export latqgp=1
export na60mode=2

export random=$(printf '%d' 0x$(xxd -l 3 -ps -c 10 /dev/urandom))
# 'relevant' parameters
iopath=HYDROPATH/EVNUM

export dt=0.2
export dxdydz=0.37

cd ${SOURCE}
export inputFile='for_dilrates_ENERGY-EVNUM.dat'
export ftn71="dil_ENERGY-EVNUM.dat"

ln ${iopath}/'for_dilrates.dat' -s -T ${inputFile}
ln ${iopath}/'dileptons_hydro.dat' -s -T ${ftn71}
singularity exec ${CON} ${SOURCE}/vhlle_read > ${iopath}/dilepton.log

rm ${inputFile}
rm ${ftn71}

exit
