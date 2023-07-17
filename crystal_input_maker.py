#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 09:56:03 2023

@author: sirous and roberto
"""
###############################################################################
##                            Modules                                       ###
###############################################################################
from pathlib import Path
import re
import numpy as np
from collections import defaultdict
from periodictable import elements
from mendeleev import element
import pickle
import time
###############################################################################
###                                                                         ###
###   Function to read lattice vector coordinates and geometries of atoms   ###
###############################################################################
def lat_vec_geom_extractor(path_to_file):
    """
    This function reads the lattice vector coordinates and geometries of atom
    from VASP

    Parameters
    ----------
    path_to_file : an Path object
        This is the path to our target file

    Returns
    -------
    3 lattice vectors as numpy array and geometries of atoms either in direct of
    cartesian 

    """
    # array to keep lattice vector and coordinates
    coord_dict = defaultdict(list)
    # pattern to read lattice vector coordinates and  xyz
    pattern = r"\s*([-+]?\d+\.\d+)\s*([-+]?\d+\.\d+)\s*([-+]?\d+\.\d+).*"
    pattern = re.compile(pattern)
    
    # we now read  our file
    with path_to_file.open() as source:
        for line in source.readlines():
            if match := pattern.match(line):
                coord_dict["X"].append(float(match.group(1)))
                coord_dict["Y"].append(float(match.group(2)))
                coord_dict["Z"].append(float(match.group(3)))
                
    # convert coord_dict to numpy array
    num_of_atom_plus_lat_vec = len(coord_dict['X'])
    # empty array to keep lattice vectors and atom coordinates
    full_array = np.zeros((num_of_atom_plus_lat_vec, 3), dtype=np.float32)
    
    #we create an array to store the lattice vectors and atoms coordinates
    #each column is given respectively to x,y,z
    for idx, key in enumerate(coord_dict.keys()):
        full_array[:,idx] = coord_dict[key]
        
    #we create separate full array into lattice vector and atoms coordinates arrays
    #lattice vector array
    latt_vect_arr = full_array[0:3]
    #atoms coordinates array
    coord_of_atoms = full_array[3:]
    
    return latt_vect_arr, coord_of_atoms

###############################################################################
###                       Atom type extractor                               ###
###############################################################################
def atom_type_extractor(path_to_file):
    """
    This function extract atom types and their numbers from our target file

    Parameters
    ----------
    path_to_file : an Path object
        This is the path to our target file

    Returns
    -------
    atom type list and their number list

    """
    
    # pattern to search
    pattern = r"\s*[^D0-9][\w|\w\w]?\s.*"
    pattern = re.compile(pattern)
    # opening and searching
    with path_to_file.open() as source:
        source_read = source.readlines()
        # ignore first line
        for line_idx in range(1, len(source_read)):
            if match := pattern.match(source_read[line_idx]):
                my_match = match.group(0).split()
                
    # empty list for storing atoms
    #atom_list = []

    
    # search on the string to find individual atoms
    #for atom_string in my_match:
    #    for let_idx in range(len(atom_string)):
    #        if let_idx != len(atom_string) - 1:
    #            
    #            if atom_string[let_idx].isupper() and atom_string[let_idx+1].isupper():
    #                atom_list.append(atom_string[let_idx])
    #            if atom_string[let_idx].isupper() and atom_string[let_idx+1].islower():
    #                atom_list.append(atom_string[let_idx]+atom_string[let_idx+1])
    #        if let_idx == len(atom_string) - 1:
    #            if atom_string[let_idx].isupper():
    #                atom_list.append(atom_string[let_idx])
                    
    return my_match
###############################################################################
###                                                                         ###
###               Extract number of each individual atoms                   ###
###                                                                         ###
###############################################################################
 # [24, 12] in integer",
def number_of_atoms_extractor(path_to_file):
     """
     This function extract number of the differernt species of atoms from our target file",

     Parameters,
     ----------,
     path_to_file : a Path object,
         This is the path to our target file,
     Returns,
     ------
     """
    #list with the number of the differernt species of atoms
     # pattern to search,
     #pattern = r"\s*[(\d+)\s*]{1,}[.]*"
     pattern = r"\s*\d+\s+\d*.*"
     pattern = re.compile(pattern)
     # opening and searching\n",
     with path_to_file.open() as source:
         for line in source.readlines():
             if match := pattern.match(line):
                 my_match = match.group()

     atom_number_list = [int(item) for item in my_match.split()]
     return atom_number_list
###############################################################################
###                                                                         ###         
###                  Atom type atomic number dictionary                     ###
###                                                                         ###
###############################################################################
periodic_table_Z_dict = defaultdict(int)

for atom in elements:
    if str(atom)!='n' :
        #print(type((element(str(atom)).econf)))
        periodic_table_Z_dict[atom.symbol] = [atom.number]
        
        
###############################################################################
###                                                                         ###         
###         Elements in electronic configuration dictionary                 ###
###                                                                         ###
###############################################################################


def electronic_configurator():
    
        """
         This function writes in a text file a dictionary in which each element 
          isassociated to a dictionary containing in turn the number of 
         electrons for each momentum-type orbital",
    
         Parameters,
         ----------,
         none
         
         Returns,
         ------
         none
         
        """
        elect_conf_dictionary=defaultdict()
        pattern = r'\s*(\d)\s(\d)\s(\d)\s(\d)\s(\d\.\d)'
        pattern_2 = r'\d(\w)(\d\d*)'
        pattern_3 = r'\d(\w)'
        pattern = re.compile(pattern)
        pattern_2 = re.compile(pattern_2)
        pattern_3 = re.compile(pattern_3)

            
        for atom in elements:
                if str(atom) != 'n' :   
            
                    period_elem_dict = {1:{'s':0,'p':0,'d':0,'f':0},
                                      2:{'s':2,'p':0,'d':0,'f':0},
                                          3:{'s':4,'p':6,'d':0,'f':0}, 
                                          4:{'s':6,'p':12,'d':0,'f':0},
                                          5:{'s':8,'p':18,'d':10,'f':0},
                                          6:{'s':10,'p':24,'d':20,'f':0},
                                          7:{'s':12,'p':30,'d':30,'f':14},
                                          8:{'s':14,'p':36,'d':40,'f':28}}
                   
                      #electron configuration retrieved form mendeelev library
                    mendel_conf = (element(str(atom)).econf).split()
                      #period of our atom
                    period = int(mendel_conf[-1][0])
                          #print(period,str(atom))
                          
                    for shell_idx in range(0,len(mendel_conf)):
                              
                              #print(elect_conf_dictionary)
                              if match := pattern_2.match(mendel_conf[shell_idx]):
                
                                  
                                  period_elem_dict[period][str(match.group(1))] += (float(
                                      match.group(2))-1)
                
                         
                              if match:=pattern_3.match(mendel_conf[shell_idx]):
                                 period_elem_dict[period][str(match.group(1))] += 1
                          
                
                    
                    elect_conf_dictionary[str(atom)]=period_elem_dict[period]
                    
        return elect_conf_dictionary
###############################################################################
###            Norms of the lattice vectors and their angles                ###
###############################################################################
def vect_norm_and_angles(latt_vect):
    """
    This function returns the norm of the lattice vectors and the angles
    between them in a single array

     Parameters,
     ----------,
     the matrix with the cartesian components of the lattice vectors
         
     Returns,
     ----------,
     A vector where the first three entries are the lattice vectors norms
     and the second one are the angles between the lattice vectors 
     (alpha,beta,gamma)
         """    
    
    #creating the vectors norms array 
    angle_arr = np.zeros(latt_vect.shape[0])
    norms_arr = np.zeros(latt_vect.shape[0])
    
    for i in range(len(latt_vect)):
        norms_arr[i]=np.linalg.norm(latt_vect[i])
        for j in range(i+1, len(latt_vect)):
                nomin = latt_vect[i]@latt_vect[j]
                domin = norms_arr[i]*np.linalg.norm(latt_vect[j])
                theta = np.arccos(nomin/domin)*180/np.pi
                angle_arr[i+j-1] = theta 

    cell_vector=np.concatenate((norms_arr,angle_arr),axis=0)
    cell_vector[3],cell_vector[5] = cell_vector[5], cell_vector[3]
    return cell_vector

###############################################################################
####                     Getting the basis set                              ###
###############################################################################
def basis_set_extractor(atom_list,basis_set_list):
    
    
    """
    This function extract the chosen basis sets for the elements extracted from
    the POSCAR file. 
    The order of the basis sets in the list given as parameters must 
    correspond to the order of the atoms in the atom species list.
    It also reads the basis set stored in the basis set dictionary
    returned by the function "basis_sets_extractor" and in the case the
    fractional charges of the electornic shells are all zeros, it assign
    to them values respect to their angular momentum.

     Parameters,
     ----------,
     1)the list with the atoms species in the POSCAR
     2)list with the several basis sets chosen for each atoms
         
     Returns,
     ----------,
    This function returns a list with the complete basis set to be inserted in
    the input
         """   
            
    #extraction of the basis set library path
    basis_set_lib_path=Path.cwd()/"BSE_library"
    #patter to recognize the shell specification
    pattern_shell = r'(\d)\s(\d)\s(\d+)\s(\d)\s(\d\.\d)'
    pattern_shell = re.compile(pattern_shell)
    #pattern to recognize the gaussians parameteres
    pattern_gauss = r'\s+\d{0,}\..+\d+'
    pattern_gauss = re.compile(pattern_gauss)
    #pattern to recognize the atom type and the number of primitives
    pattern_prim = r'(\d+\s\d+)[\n]'
    pattern_prim = re.compile(pattern_prim)
    #creation of the storing list
    basis_list=[]
    #dictionary to convert each shell type in the corresponding CRYSTAL number
    conversion_dictionary = {'s':0,'p':2,'d':3,'sp':1,}
    #we look for the atom in the atoms list
    for i in range(len(atom_list)):
        #temporary dictionary to store all electronic configuration for each atom
        temp_dict=defaultdict()
        #reasearch of the atom folder in the library
        for file in basis_set_lib_path.iterdir():
            if atom_list[i] == file.name:
                #filling of the temporary dictionary with the electrons brought
                #the electronic configuration database
                for key in conversion_dictionary:
                    #assignment for shell s and p
                    if key != 'd' and key !='sp':
                        temp_dict[conversion_dictionary[key]] = electronic_config_database_loaded[atom_list[i]][key]
                    #assignment for shell d
                    elif key == 'd':
                        temp_dict[conversion_dictionary[key]] = electronic_config_database_loaded[atom_list[i]][key]
                        +electronic_config_database_loaded[atom_list[i]]['f']
                    #assigment for shell sp
                    else :
                        temp_dict[conversion_dictionary[key]] = electronic_config_database_loaded[atom_list[i]]['s']+ electronic_config_database_loaded[atom_list[i]]['p'] 
                            
                #research of the specified basis set in the atom folder
                for j in file.iterdir():
                    #stem attribute takes the last part of the path soecified in j
                    if j.stem == basis_set_list[i]:
                        with j.open() as source:
                                #starting to read each line of the basis set
                                for line in source.readlines():
                                   if match := pattern_shell.match(line):
                                       shell_type = int(match.group(2))
                                       shell_num = int(match.group(4))
                                       #exclusion of shells f and g and check of there are available electrons to put
                                       if  not shell_type in [4,5,6,7] and temp_dict[shell_type] > 0 :
                                           if shell_type == 0 and temp_dict[shell_type] > 1 or shell_type == 1 and temp_dict[shell_type] > 1:
                                               shell_num += 2 
                                               temp_dict[0] -= shell_num
                                               temp_dict[1] -= shell_num
                                           elif shell_type == 0 or shell_type == 1 :
                                               shell_num += temp_dict[shell_type]
                                               temp_dict[0] -= shell_num
                                               temp_dict[1] -= shell_num
                                           if shell_type == 2 and temp_dict[shell_type] >= 6 or shell_type == 1 and temp_dict[shell_type] >= 6:
                                               shell_num += 6
                                               temp_dict[2] -= shell_num
                                               temp_dict[1] -= shell_num
                                           elif shell_type == 2 or shell_type == 1 :
                                               shell_num += temp_dict[shell_type]
                                               temp_dict[2] -= shell_num
                                               temp_dict[1] -= shell_num
                                           if shell_type == 3 and temp_dict[shell_type] >= 10:
                                               shell_num += 10
                                               temp_dict[3] -= shell_num
                                           elif shell_type == 3: 
                                               shell_num += temp_dict[shell_type]
                                               temp_dict[3] -= shell_num    
                                       #string with all the matches from the pattern_shell with the upgraded number of electrons
                                       #per shell, it fills the shell_list
                                       basis_list.append(f"{match.group(1)} {match.group(2)} {match.group(3)} {shell_num} {match.group(5)}")
                                   if match := pattern_gauss.match(line) :
                                       basis_list.append(match.group(0))
                                   if match := pattern_prim.match(line) :
                                       basis_list.append(match.group(1))
                        break
                print(basis_list)
                
                break         
    basis_list.remove('99 0') 
    basis_list.remove('99 0')
    basis_list.append('99 0')


    
    return basis_list 
###############################################################################
###                                                                         ###
###                         Call electronic_configurator                    ###
###                                                                         ### 
###############################################################################
#electronic_configuration_dict = electronic_configurator()
# save electronic configuration database to a permanent file
t1_1 = time.time()
#with open('electronic_config_database.pkl', 'wb') as source:
#    pickle.dump(electronic_configuration_dict, source)
#    print('dictionary saved successfully to file')
#t1_2 = time.time()
#print(f"writing time: {t1_2-t1_1}")   
# read the dictionary from our save dataset
#t2_1 = time.time()
with open('electronic_config_database.pkl', 'rb') as source:
    electronic_config_database_loaded = pickle.load(source)
#t2_2 = time.time()
#print(f"reading time: {t2_2-t1_2}")   
#print(electronic_config_database_loaded)


#addition of the atomic numbers to the coordinates
###############################################################################
###                   Path to Poscar geometry                               ###
###############################################################################
# current working directory
CWD = Path.cwd()
for file in Path(CWD/"POSCAR").iterdir():
    if file.name == "IIPOSCAR.txt":
        path_to_xyz_poscar = file

latt_vect_coordinates, coord_of_atoms = lat_vec_geom_extractor(path_to_xyz_poscar)

atom_type = atom_type_extractor(path_to_xyz_poscar)

number_of_atoms = number_of_atoms_extractor(path_to_xyz_poscar)

cell_geom = vect_norm_and_angles(latt_vect_coordinates)

#electronic_configurator()



#print(t1_1, t1_2, t1_3, t1_4, t1_5, t1_6, t1_7)
###############################################################################
###                     Atomic number of our atoms!                         ###
###                                                                         ### 
###############################################################################
# empty list to store symbol of atoms
symbol_of_all_atoms_list = [] 
# empty list for atomic number of all atoms
atomic_number_list = []

for atom_idx in range(len(atom_type)):
    for _ in range(number_of_atoms[atom_idx]):
        symbol_of_all_atoms_list.append(atom_type[atom_idx])
        

for atom in symbol_of_all_atoms_list:
    atomic_number_list.append(int(periodic_table_Z_dict[atom][0]))

full_matrix_coord_crystal=np.column_stack([np.array(atomic_number_list),
                                           coord_of_atoms])

###############################################################################
###                 writing input geomtry blocks                            ###
###############################################################################

with open("_geometry_block.txt",'w') as source:
    source.write("CRYSTAL\n")
    source.write("0 0 0\n")
    source.write("1\n")
    for number in cell_geom:
        source.write(str(number)+" ")
    source.write("\n")
    source.write(str(np.array(number_of_atoms).sum())+"\n")
    for coordinates in full_matrix_coord_crystal:
            for numbers in list(coordinates):
                source.write(f"{numbers} ")
            source.write("\n")
    
###############################################################################
###                     writing input basis set blocks                      ###
###############################################################################

sets=[]
with Path.cwd().joinpath("basis_block_new","sets_list.txt") as source:
     sets=re.split("\s+|\n",source.read_text())[:-1]
    
print(sets)
  
for i in sets:

     dir = Path.cwd().joinpath("basis_block_new",i)

     try:
        base= basis_set_extractor(atom_type,[i,i])
        with dir.with_suffix(".txt").open(mode='w') as source:
            for line in base:
                source.write(str(line)+'\n')
     except:
        print(' this basis set is not working', i)
    
    
    


    
    
    
             


        