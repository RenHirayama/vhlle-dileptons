# Fortran compiler - uncomment the one you want

 F77 = gfortran -ffixed-line-length-132
# F77 = ifort -132

# Flags for Fortran Compilers:  Pick one for each compiler and uncomment it.

#### FOR gfortran ####
# F77FLAGS = -O2 -mcmodel=medium -falign-commons -fforce-addr -funroll-loops #-fbacktrace
 F77FLAGS = -O3 -mcmodel=medium #-falign-commons -fforce-addr -funroll-loops #-fbacktrace
#### FOR ifort ####
# F77FLAGS = -O3 -ipo -fcommon -mcmodel=medium


#coarse			:coarse.o input.o readeos.o grid.o transform.o temp.o freezeout.o 4piemit.o 4piemit_mix.o dilemit.o dilemit_rapp.o dilemit_rapp_hr.o qgpemit.o qgpemit_lat.o photon.o dirrho.o omega.o dirphi.o shdelta.o pel_ppo_dist.o output.o ranf.o pdg2ityp.o
#			$(F77) $(F77FLAGS) -o coarse coarse.o input.o readeos.o grid.o transform.o temp.o freezeout.o 4piemit.o 4piemit_mix.o dilemit.o dilemit_rapp.o dilemit_rapp_hr.o qgpemit.o qgpemit_lat.o photon.o dirrho.o omega.o dirphi.o shdelta.o pel_ppo_dist.o output.o ranf.o pdg2ityp.o
#
#coarse.o		:coarse.f
#			$(F77) $(F77FLAGS) -c coarse.f
#input.o			:input.f 
#			$(F77) $(F77FLAGS) -c input.f

vhlle_dil			:vhlle_dil.o vhlle_read.o readeos.o grid.o transform.o temp.o freezeout.o 4piemit.o 4piemit_mix.o dilemit.o dilemit_rapp.o dilemit_rapp_hr.o qgpemit.o qgpemit_lat.o photon.o dirrho.o omega.o dirphi.o shdelta.o pel_ppo_dist.o output.o ranf.o pdg2ityp.o
			$(F77) $(F77FLAGS) -o vhlle_dil vhlle_dil.o vhlle_read.o readeos.o grid.o transform.o temp.o freezeout.o 4piemit.o 4piemit_mix.o dilemit.o dilemit_rapp.o dilemit_rapp_hr.o qgpemit.o qgpemit_lat.o photon.o dirrho.o omega.o dirphi.o shdelta.o pel_ppo_dist.o output.o ranf.o pdg2ityp.o

vhlle_dil.o		:vhlle_dil.f
			$(F77) $(F77FLAGS) -c vhlle_dil.f
vhlle_read.o			:vhlle_read.f 
			$(F77) $(F77FLAGS) -c vhlle_read.f
readeos.o		:readeos.f 
			$(F77) $(F77FLAGS) -c readeos.f
grid.o			:grid.f
			$(F77) $(F77FLAGS) -c grid.f
transform.o		:transform.f
			$(F77) $(F77FLAGS) -c transform.f
temp.o			:temp.f
			$(F77) $(F77FLAGS) -c temp.f
freezeout.o		:freezeout.f
			$(F77) $(F77FLAGS) -c freezeout.f
4piemit.o		:4piemit.f
			$(F77) $(F77FLAGS) -c 4piemit.f
4piemit_mix.o		:4piemit_mix.f
			$(F77) $(F77FLAGS) -c 4piemit_mix.f
dilemit.o		:dilemit.f
			$(F77) $(F77FLAGS) -c dilemit.f
dilemit_rapp.o		:dilemit_rapp.f
			$(F77) $(F77FLAGS) -c dilemit_rapp.f
dilemit_rapp_hr.o	:dilemit_rapp_hr.f
			$(F77) $(F77FLAGS) -c dilemit_rapp_hr.f
qgpemit.o		:qgpemit.f
			$(F77) $(F77FLAGS) -c qgpemit.f
qgpemit_lat.o		:qgpemit_lat.f
			$(F77) $(F77FLAGS) -c qgpemit_lat.f
photon.o		:photon.f
			$(F77) $(F77FLAGS) -c photon.f
dirrho.o		:dirrho.f
			$(F77) $(F77FLAGS) -c dirrho.f
omega.o			:omega.f
			$(F77) $(F77FLAGS) -c omega.f
dirphi.o		:dirphi.f
			$(F77) $(F77FLAGS) -c dirphi.f
shdelta.o		:shdelta.f
			$(F77) $(F77FLAGS) -c shdelta.f
pel_ppo_dist.o		:pel_ppo_dist.f
			$(F77) $(F77FLAGS) -c pel_ppo_dist.f
output.o		:output.f
			$(F77) $(F77FLAGS) -c output.f
ranf.o			:ranf.f
			$(F77) $(F77FLAGS) -c ranf.f
pdg2ityp.o		:pdg2ityp.f
			$(F77) $(F77FLAGS) -c pdg2ityp.f

clear			:			
			rm *.o
