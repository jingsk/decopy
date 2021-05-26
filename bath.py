#!/usr/bin/env python
# coding: utf-8

# In[5]:


import numpy as np
import pandas as pd
from ase.io import read
from ase.visualize import view
from wrapcube import cleanUp
from ase.data.isotopes import download_isotope_data
#use like this isotopes[6][12]['composition'] gives 0.989
isotopes = download_isotope_data()
from ase.symbols import string2symbols, symbols2numbers


#usage guide
#SWCNT=bath('CONTCAR')
#SWCNT.bath_geometry

class bath:
    def __init__(self, fname: str,supercell: tuple = (1,1,1), r_cutoff: float=0.0):
        #wrap and center if vacuum exists
        #TODO: create a clean supercell with an origin at the center (low priority)
        self.atoms=cleanUp(read(fname))*supercell
        #for atom species in atoms get isotope abundance data
        self.spin_table = get_spin_table(self.atoms)
        #randomly generate spin bath geometry based on % abundance
        #df of symbol (~p), x, y, z, isotope (n+p), nuclear_spin (I), distance
        self.bath_geometry= generate_spin_sites(self.atoms,self.spin_table)
        #TODO: determine spin site (sp3 carbon site 
        #(is this possible or do we need to ask user for e- spin site?))
        
#for present atomic species get isotope data and return
#a DataFrame with atomic number, atomic mass, percent abundance, nuclear_spin
def get_spin_table(atoms):
    isotope_dfs=[]
    for atom in get_unique_atomic_species(atoms):
        atom_isotope_df=pd.DataFrame(isotopes[atom]).T
        # add atomic_number
        atom_isotope_df.loc[:,'atomic_number']=atom
        #remove unnatural isotopes 
        atom_isotope_df=atom_isotope_df[atom_isotope_df.composition>0]
        #I=(n-p)/2
        atom_isotope_df.loc[:,'nuclear_spin']=get_nuclear_spin(atom_isotope_df.index-atom,atom)
        isotope_dfs.append(atom_isotope_df)
    return pd.concat(isotope_dfs)

#return a DataFrame of spin active sites inside a cutoff radius
#TODO: add an option to input a cutoff radius, add input for specifying a spin site as a point
#DataFrame columns include: atomic symbol, isotope,
#x,y,z, nuclear spin(I), and distance from a point
def generate_spin_sites(atoms,spin_table):
    #make xyz-format df
    positions=pd.concat([pd.DataFrame(atoms.get_chemical_symbols(),columns=['symbol']), 
               pd.DataFrame(atoms.get_positions(),columns=['x','y','z'])],
               axis=1)
    spin_active_positions_list=[]
    for species in np.unique(positions.symbol):
        #get atomic number
        atomic_number=symbols2numbers(species)[0]
        #feed in spin_table of the atomic species
        #and dataframe of positions of the atomic species
        #to assign randomly generated isotope
        species_positions=assign_isotopes(positions[positions.symbol==species],spin_table[spin_table.atomic_number==atomic_number])
        #calculate nuclear spin from (p-n)/2
        species_positions.loc[:,'nuclear_spin']=get_nuclear_spin(species_positions.isotope-atomic_number,atomic_number)
        #remove nuclear spin inactive sites
        species_positions=species_positions[is_spin_active(species_positions.nuclear_spin)]
        #find distance from electron spin site (assuming at 0,0,0)
        species_positions.loc[:,'distance']=get_distance_from_point(species_positions[['x','y','z']])
        #remove sites beyond r_cutoff
        species_positions=species_positions[within_r_cutoff(species_positions.distance)]
        spin_active_positions_list.append(species_positions)
    return pd.concat(spin_active_positions_list)


def get_unique_atomic_species(atoms):
    return np.unique(atoms.get_atomic_numbers())
def get_nuclear_spin(num_p,num_n):
    return np.abs((num_n-num_p)/2)
def is_spin_active(nuclear_spin):
    return nuclear_spin!=0
def get_distance_from_point(positions,point=[0,0,0]):
    rel_x=positions.x-point[0]
    rel_y=positions.y-point[1]
    rel_z=positions.z-point[2]
    return np.sqrt(rel_x**2+rel_y**2+rel_z**2)
#r_cutoff in angstrom
def within_r_cutoff(distance,r_cutoff=100):
    return distance <r_cutoff

def assign_isotopes(positions,isotope_df):
    #create randomly generated numbers from 0 to 1
    positions.loc[:,'rand']=np.random.rand(positions.shape[0])
    #generate isotope condition
    #this works like this: if X-12, X-13, X-14 has [0.1,0.2,0.7] composition
    #starting with all X-12, we assign isotope with rand>0.1 to X-13, rand>0.3 to X-14 
    isotope_condition=np.cumsum(np.array(isotope_df.composition))
    #initial assignment
    positions.loc[:,'isotope']=isotope_df.index[0]
    #in case there's only one isotope then do nothing
    if isotope_df.index.size >1:
        #assign isotope based on abundance condition
        for i in np.arange(1,isotope_df.index.size):
            positions.loc[positions.rand>isotope_condition[i-1], "isotope"]=isotope_df.index[i]
    positions.drop('rand', axis=1, inplace=True)
    return positions


SWCNT=bath('CONTCAR_NO2',(1,1,3))
SWCNT.bath_geometry


# In[2]:


#view(SWCNT.atoms)


# In[3]:


#atoms=read('CONTCAR_NO2')
#atoms*(1,1,3)


# In[4]:


#atoms.positions

