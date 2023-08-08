# SOP-HYBRID
SOP-HYBRID is a LAMMPS package for simulating IDP and Multidomain Proteins


Usage instructions:
Place this folder in the src folder of the LAMMPS source code and use
$ make yes-SOP-HYBRID

## In Serial environments
$ make serial
## In Parallel environments
$ make mpi

Once successfully compiled you will notice the executable lmp_serial or lmp_mpi. This executable can be used along with the input files
to run the simulations.

This Package has been tested with lammps-23Jun2022 version
