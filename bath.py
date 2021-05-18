import numpy as np
import pandas as pd
from ase.io import read
from ase.visualize import view
from wrapcube import cleanUp
from ase.data.isotopes import download_isotope_data
from ase.symbols import string2symbols, symbols2numbers
isotopes = download_isotope_data()

#usage guide
#SWCNT=bath('CONTCAR')
#SWCNT.bath_geometry

class bath:
    #fname is file name to read
    #supercell a tuple for repeating unit cell in x,y,z 
    def __init__(self, fname,supercell, r_cutoff=0.0):
        #wrap and center if vacuum exists
        self.atoms=cleanUp(read(fname))*supercell
        #for atom species in atoms get isotope abundance data
        self.spin_table = get_spin_table(self.atoms)
        #randomly generate spin bath geometry based on % abundance
        #df of symbol (~p), x, y, z, isotope (n+p), nuclear_spin (I), distance
        self.bath_geometry= generate_spin_sites(self.atoms,self.spin_table)
        
    #
def get_spin_table(atoms):
#isotope_df=pd.DataFrame({})
    #atoms=self.atoms
    isotope_dfs=[]
    for atom in get_unique_atomic_species(atoms):
        atom_isotope_df=pd.DataFrame(isotopes[atom]).T
        atom_isotope_df['atomic_number']=atom
        atom_isotope_df=atom_isotope_df[atom_isotope_df.composition>0]
        #isotope_df.index.name='num_p+n'
        atom_isotope_df['nuclear_spin']=get_nuclear_spin(atom_isotope_df.index-atom,atom)
        #atom_isotope_df=atom_isotope_df[is_spin_active(atom_isotope_df.nuclear_spin)]
        isotope_dfs.append(atom_isotope_df)
    return pd.concat(isotope_dfs)
def generate_spin_sites(atoms,spin_table):
    #atoms=self.atoms
    #spin_table=self.spin_table
    #make xyz like df
    positions=pd.concat([pd.DataFrame(atoms.get_chemical_symbols(),columns=['symbol']), 
               pd.DataFrame(atoms.get_positions(),columns=['x','y','z'])],
               axis=1)
    spin_active_positions_list=[]
    #loop through elements in atomic positions
    for species in np.unique(positions.symbol):
        #dataframe of positions of the atomic species
        species_positions=positions[positions.symbol==species]
        #get atomic number
        atomic_number=symbols2numbers(species)[0]
        #spin_table of the atomic species
        atom_isotope_df=spin_table[spin_table.atomic_number==atomic_number]
        #create a column in the df of randomly generated numbers from 0 to 1
        species_positions['rand']=np.random.rand(species_positions.shape[0])
        #generate isotope condition
        #this works like this: if X-12, X-13, X-14 has [0.1,0.2,0.7] composition
        #starting with all X-12, we assign isotope with rand>0.1 to X-13, rand>0.3 to X-14 
        isotope_condition=np.cumsum(np.array(atom_isotope_df.composition))
        #initial assignment
        species_positions['isotope']=atom_isotope_df.index[0]
        #in case there's only one isotope then do nothing
        if atom_isotope_df.index.size >1:
            #assign isotope based on abundance condition
            for i in np.arange(1,atom_isotope_df.index.size):
                species_positions.loc[species_positions.rand>isotope_condition[i-1], "isotope"]=atom_isotope_df.index[i]
        species_positions.drop('rand', axis=1, inplace=True)
        #calculate nuclear spin from (p-n)/2
        species_positions['nuclear_spin']=get_nuclear_spin(species_positions.isotope-atomic_number,atomic_number)
        #remove nuclear spin inactive sites
        species_positions=species_positions[is_spin_active(species_positions.nuclear_spin)]
        #find distance from electron spin site (assuming at 0,0,0)
        species_positions['distance']=get_distance_from_point(species_positions[['x','y','z']])
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