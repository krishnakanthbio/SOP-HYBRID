#	To generate a .pdb file with the provided sequence
#	Krishnakanth B, 
#	Theoretical Biophysics Laboratory, Molecular Biophysics Unit,
#	Indian Institute of Science, Bangalore - 560012
#
#	Last Modified: 14 March 2023

#  Load the required packages

import __main__
__main__.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI

import sys, time, os, cmd
import pymol

pymol.finish_launching()

##
# Read User Input
#spath = os.path.abspath(sys.argv[1])
#sname = spath.split('/')[-1].split('.')[0]

# Load Structures

#pymol.cmd.load(spath, sname)
#pymol.cmd.disable("all")
#pymol.cmd.enable(sname)
#pymol.cmd.png("my_image.png")

# Read user sequence from provided fasta file
try:
    sys.argv[1]
except IndexError as ie:
    print("Exception : {0}".format(ie))


fasta_file = open(sys.argv[1],'r').readlines()[1:]
fasta_seq = ""

for i in fasta_file:
	fasta_seq = fasta_seq + i.strip()
	if i.strip() == "":
		break

for aa in fasta_seq: pymol.cmd._alt(str.lower(aa))

pdb_fname = sys.argv[1].split(".fasta")[0]+".pdb"
pymol.cmd.save(pdb_fname, selection='(all)', state=-1, format='', ref='',
         ref_state=-1, quiet=1, partial=0)



# Get out!
pymol.cmd.quit()

