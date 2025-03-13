      real*8 pi,pi2,hc3,e_0,n_0,Bag,alpha_em,hqc
      real*8 dpt,dm,dmu,dptr,dmr,dmr_lr,dmr_hr,dmur,dtempr,dqr
      real*8 mass_electron,mass_muon,mass_pion,m_lept,mmin,mmax
      real*8 rhoeffmax
      integer imumax,jtmax,kqmax,kmmax
      integer itmax,jrhomax,kpimax,lkmax,nqmax,mmmax
      integer itmaxor,jrhomaxor,kpimaxor,lkmaxor,nqmaxor,mmmaxor
      integer itmaxrhom,jrhomaxrhom,kpimaxrhom,lkmaxrhom,
     &        nqmaxrhom,mmmaxrhom
      integer itmaxphi,jrhomaxphi,kpimaxphi,lkmaxphi,nqmaxphi,mmmaxphi
      integer evtmax,outmax,timemax

      integer eos
      integer seedinit
      integer rates
      integer pikchem
      integer baryons
      integer fourpimix
      integer outputs
      logical dimuon
      logical ext_out
      logical htlcorr
      logical focalc
      integer latqgp
      integer na60mode
      integer hgpar
      integer phqgp
      integer model

      logical vHLLE_out

c-----|----------------------------------------------------------------|X
c. Define constants

      PARAMETER (pi=3.141592653589793d0)
      PARAMETER (pi2=(3.141592653589793d0)**2)
      PARAMETER (hc3=197.327d0**3) ! value of hbar*c **3 in (Mev*fm)**3
      PARAMETER (e_0=146.51751415742d0)
      PARAMETER (n_0=0.15891d0)
      PARAMETER (Bag=235d0**4)
      PARAMETER (alpha_em=1.d0/137.d0)
      PARAMETER (mass_electron=0.000511d0)
      PARAMETER (mass_muon=0.1057d0)
      PARAMETER (mass_pion=0.139d0)
      PARAMETER (hqc=0.197327d0) ! value of hbar*c in GeV*fm

c-----|----------------------------------------------------------------|X
c. ARRAY SIZES
c.. LHC Energies:
c      parameter (evtmax=1000)
c      parameter (timemax=126)
c      parameter (outmax=50000)
c.. RHIC Energies:
c      parameter (evtmax=1000)
c      parameter (timemax=126)
c      parameter (outmax=16500)
c.. NA60 Energies:
c      parameter (evtmax=1001)
c      parameter (timemax=126)
c      parameter (outmax=10000)
c.. SIS Energies
      parameter (evtmax=1600)
      parameter (timemax=600)
      parameter (outmax=500)

c-----|----------------------------------------------------------------|X
c. For dilepton table read-in:

c.. Bins of the dilepton table (ELETSKY)
c..   PARAMETER (imumax=200,jtmax=301,kqmax=300,kmmax=150)
c..   PARAMETER (imumax=200,jtmax=201,kqmax=300,kmmax=150)
      PARAMETER (imumax=200,jtmax=271,kqmax=300,kmmax=150)

c.. Bins of the dilepton table (RAPP)
      PARAMETER (itmax=32,jrhomax=28,kpimax=4,lkmax=4)
      PARAMETER (nqmax=21,mmmax=31)

c.. Bins of the dilepton table (RAPP - RHO & OMEGA ONLY)
      PARAMETER (itmaxrhom=32,jrhomaxrhom=13,kpimaxrhom=4,lkmaxrhom=1)
      PARAMETER (nqmaxrhom=21,mmmaxrhom=150)

c.. Bins of the dilepton table (RAPP - PHI ONLY)
      PARAMETER (itmaxphi=32,jrhomaxphi=13,kpimaxphi=4,lkmaxphi=4)
      PARAMETER (nqmaxphi=21,mmmaxphi=101)

c. ELETSKY RESOLUTIONS:
c.. Resolution of the pt-grid
      PARAMETER (dpt=5.d-3) !5 MeV
c.. Resolution of the mass-grid
      PARAMETER (dm=1.d-2) !10 MeV
c.. Resolution of th mub-grid
      PARAMETER (dmu=1.d-2) !10 MeV

c. RAPP RESOLUTIONS:
c.. Resolution of the momentum-grid
      PARAMETER (dqr=25.d-2) !250 MeV
c.. Resolution of the pt-grid
      PARAMETER (dptr=25.d-2) !255 MeV
c.. Resolution of the mass-grid
      PARAMETER (dmr_lr=5.d-2) !50 MeV
      PARAMETER (dmr_hr=1.d-2) !50 MeV
c.. Resolution of the rho-grid
      PARAMETER (dmur=5.d-1) !0.5 rho_0
c.. Resolution of the temp-grid
      PARAMETER (dtempr=1.d-2) !10 MeV

c-----|----------------------------------------------------------------|X
c.  ***Choose***:
c.. Output format
      PARAMETER(vHLLE_out=.TRUE.) ! Standard output
      PARAMETER(ext_out=.FALSE.) ! Standard output
c      PARAMETER(ext_out=.TRUE.) ! Extended output

c-----|----------------------------------------------------------------|X
c.  ***Choose***:
c.. Hard-thermal loop (HTL) correction for QGP rates (does not apply
c.. for Lattice Rates)
c      PARAMETER(htlcorr=.FALSE.) ! No HTL correction (Pure QED)
      PARAMETER(htlcorr=.TRUE.) ! With HTL correction

c-----|----------------------------------------------------------------|X

c.. Minimum/Maximum mass
        common/mthres/ mmin,mmax
c..Bins of the dilepton table
        common/massres/ dmr,rhoeffmax
        common/tableres/ itmaxor,jrhomaxor,kpimaxor,lkmaxor,nqmaxor,
     &                   mmmaxor

c-----|----------------------------------------------------------------|X
        common /randomf/ seedinit
        common /inputs/  eos,dimuon,rates,pikchem,baryons,
     &                   fourpimix,outputs,latqgp,na60mode,focalc,model
        common /inputph/ hgpar,phqgp
