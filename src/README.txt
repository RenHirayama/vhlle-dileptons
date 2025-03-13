# Coarse Graining by Stephan Endres

## to build project

run:
> make

## to run

* configure by modifying *runfile*
* input and output also defined here
* IMPORTANT for SMASH: define e_cm and beta_lab

Werte für Gold+Gold bei 1.23 AGeV sind

E_lab(GeV/u): 0.1230E+01
sqrt(s)(GeV): 0.2414E+01
beta lab fr : 0.6292914

* SMASH output that is used: particle_lists.oscar with intermediate output (e.g. every 0.7 fm)

run:
> ./runfile

## to run on cluster

use:
> sbatch ./submit_runfile

## for analysis

* see analysis folder
* also output_formats.txt


(see also Notebook I, p.97 + Notebook II, p.25)
