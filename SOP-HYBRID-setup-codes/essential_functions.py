#	Essential supporting functions required in other codes
#	Krishnakanth B, 
#	Theoretical Biophysics Laboratory, Molecular Biophysics Unit,
#	Indian Institute of Science, Bangalore - 560012
#
#	Last Modified:  28 April 2023

import numpy as np

# Euclidean distance rounded to 4 digits after decimal
def euclidean_distance(a, b):
	return round(((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)**0.5,4)



#Function to read a .pdb file, obeys 1996 nomenclature of .pdb file
def pdb_reader(pdb_filename):
	print("Reading pdb file: "+pdb_filename)
	x, y, z,res_id,atom_id ,atom_type,res_type,res_no,chain_id = [],[],[],[],[],[],[],[],[]
	f=open(pdb_filename,"r").readlines()
	for i in range(len(f)):
		if(f[i][0:4]=="ATOM"):
			x.append(float(f[i][30:38]))
			y.append(float(f[i][38:46]))
			z.append(float(f[i][46:54]))
			atom_id.append(str.strip(f[i][6:11]))
			atom_type.append(str.strip(f[i][12:16]))
			res_type.append(str.strip(f[i][17:20]))
			res_no.append(str.strip(f[i][22:26]))
			chain_id.append(f[i][21])
	return (x,y,z, atom_id, atom_type, res_type, res_no, chain_id)

# Get the unique elements of a list : Used to get the total number of chains
def unique_list(my_list):
	uniq_list=[]
	for i in  my_list:
		if not (i  in uniq_list):
			uniq_list.append(i)
	return uniq_list

# which like function (R language) in python
def which(x,my_list):
	req=[]
	for i in range(len(my_list)):
		if(x==my_list[i]):
			req.append(i)
	return req
	
# This function takes in an np array of atom types and cleans the H atom types starting with numeric value
def clean_H(y):
	for i in range(len(y)):
		if(str.isnumeric(y[i][0])):
			y[i] = y[i][1:]+y[i][0]
#	return(y)
	
# Function to read a large file in chunks
def rows(f, chunksize=1024, sep='|'):
    """
    Read a file where the row separator is '|' lazily.

    Usage:

    >>> with open('big.csv') as f:
    >>>     for r in rows(f):
    >>>         process(r)
    """
    row = ''
    while (chunk := f.read(chunksize)) != '':   # End of file
        while (i := chunk.find(sep)) != -1:     # No separator found
            yield row + chunk[:i]
            chunk = chunk[i+1:]
            row = ''
        row += chunk
    yield row
	
	
def rg_ree_calc(x,y,z, atom_type, mass_dict):
	rg = 0.0
	cmx = 0.0
	cmy = 0.0
	cmz = 0.0
	t_mass = 0.0
	
	for i in range(len(x)):
		cmx = cmx + x[i]*mass_dict[atom_type[i]]
		cmy = cmy + y[i]*mass_dict[atom_type[i]]
		cmz = cmz + z[i]*mass_dict[atom_type[i]]
		t_mass = t_mass + mass_dict[atom_type[i]]
	cmx = cmx/t_mass 
	cmy = cmy/t_mass
	cmz = cmz/t_mass
	for i in range(len(x)):
		rg = rg + mass_dict[atom_type[i]]*((x[i]-cmx)**2 + (y[i]-cmy)**2 + (z[i]-cmz)**2 )
	rg = np.sqrt(rg/t_mass)
	ree = np.sqrt(((x[0]-x[-1])**2) + ((y[0]-y[-1])**2) + ((z[0]-z[-1])**2))
	return([rg, ree])
