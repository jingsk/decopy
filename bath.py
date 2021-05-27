#!/usr/bin/env python
# coding: utf-8

# In[64]:


import numpy as np
import pandas as pd
from ase import Atoms
from ase.io import read
from ase.visualize import view
from wrapcube import cleanUp
from ase.data.isotopes import download_isotope_data
#use like this isotopes[6][12]['composition'] gives 0.989
isotopes = download_isotope_data()
from ase.symbols import string2symbols, symbols2numbers

class bath:
    """
    A class to represent a spin-active nuclear bath.
    Takes in positions and randomly assigns isotopes.

    Attributes
    ----------
    atoms: 
        ase's Atoms object to get geometry
    spin_table: 
        contains atomic number, atomic mass, percent abundance, nuclear_spin
        of species present in self.atoms
    bath_geometry: 
        contains atomic symbol, isotope, xyz positions, nuclear spin(I),
        and the distance from e- Spin

    Methods
    -------
    apply_rcutoff(e_spin_position, cutoff_radius):
        remove sites beyond a certain radius from e- spin
        use when e- spin has been defined
    
    Usage Guide
    -------
        SWCNT=bath.from_file('CONTCAR_NO2',(1,1,3))
        SWCNT.bath_geometry
        SWCNT.apply_r_cutoff(origin=[0,0,0], r=40)
        SWCNT.bath_geometry
    """
    def __init__(self, atoms: ase.Atoms):
        """
        Constructs all the necessary attributes for the bath object.

        Parameters
        ----------
            atoms:
                atoms to make bath
        """
        self.atoms=atoms
        #for atom species in atoms get isotope abundance data
        self.spin_table = get_spin_table(self.atoms)
        #randomly generate spin bath geometry based on % abundance
        #df of symbol (~p), x, y, z, isotope (n+p), nuclear_spin (I), distance
        self.bath_geometry= generate_spin_sites(self.atoms,self.spin_table)
        #TODO: determine spin site (sp3 carbon site 
        #(is this possible or do we need to ask user for e- spin site?))
        
    @classmethod
    #use cleanUp function to center atoms around surrounding vacuum
    def from_file(cls, filename: str, repeat: tuple = (1, 1, 1)):
        """
        read Atoms from file (assuming cell, positions present)
        and make supercell.

        Parameters
        ----------
            filename: 
                file name
            repeat: 
                multiple to multiply self.atoms to make supercell
        """
        return cls(cleanUp(read(filename))*repeat)
    
    def apply_r_cutoff(self, origin: np.ndarray, r: float):
        """
        remove sites beyond a certain radius from e- spin
        use when e- spin has been defined

        Parameters
        ----------
            origin:
                cart. position of e- spin 
            r:
                cutoff radius
        """
        #find distance from electron spin site (assuming at 0,0,0)
        self.bath_geometry.loc[:,'distance']=get_distance_from_point(self.bath_geometry[['x','y','z']].values)
        print(within_r_cutoff(self.bath_geometry.distance, r_cutoff=r))
        #remove sites beyond r_cutoff
        self.bath_geometry=self.bath_geometry[within_r_cutoff(self.bath_geometry.distance, r_cutoff=r)].copy()

def get_spin_table(atoms: ase.Atoms):
    """
    for present atomic species get isotope data and return
    a DataFrame with atomic number, atomic mass, percent abundance, nuclear_spin
    
    Parameters
    ----------
        atoms:
            atoms with elements to search for isotopes 
    """
    isotope_dfs=[]
    #atom is atomic number
    for atomic_number in get_unique_atomic_number(atoms):
        #obtain isotope data from ase.data.isotopes
        atom_isotope_df=pd.DataFrame(isotopes[atomic_number]).T
        atom_isotope_df.loc[:,'atomic_number']=atomic_number
        #remove unnatural isotopes 
        atom_isotope_df=atom_isotope_df[atom_isotope_df.composition>0]
        #I=(n-p)/2
        atom_isotope_df.loc[:,'nuclear_spin']=get_nuclear_spin(atom_isotope_df.index-atomic_number,atomic_number)
        isotope_dfs.append(atom_isotope_df)
    return pd.concat(isotope_dfs)

#return a DataFrame of spin active sites

def generate_spin_sites(atoms: ase.Atoms,spin_table: pd.DataFrame):
    """
    gives spin active sites as DataFrame based on isotopes 
    present based on positions and atomic species in atoms
    
    Parameters
    ----------
        atoms:
            atoms with elements to search for isotopes 
    """
    #make xyz-format df
    positions=pd.concat([pd.DataFrame(atoms.get_chemical_symbols(),columns=['symbol']), 
               pd.DataFrame(atoms.get_positions(),columns=['x','y','z'])],
               axis=1)
    spin_active_positions_list=[]
    for species in np.unique(positions.symbol):
        #get atomic number of one atom
        atomic_number=symbols2numbers(species)[0]
        species_positions=assign_isotopes(positions[positions.symbol==species],spin_table[spin_table.atomic_number==atomic_number])
        #calculate nuclear spin from (p-n)/2
        species_positions.loc[:,'nuclear_spin']=get_nuclear_spin(species_positions.isotope-atomic_number,atomic_number)
        #remove nuclear spin inactive sites
        species_positions=species_positions[is_spin_active(species_positions.nuclear_spin)]
        spin_active_positions_list.append(species_positions)
    return pd.concat(spin_active_positions_list)

def assign_isotopes(positions: pd.DataFrame,isotope_df: pd.DataFrame):
    """
        feed in spin_table of the atomic species
        and dataframe of positions of the atomic species
        to assign randomly generated isotope
    
    Parameters
    ----------
        positions: 
            xyz format positions of one atomic species
        isotope_df: 
            spin_table of present isotopes of one atomic species
        
    """
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


def get_unique_atomic_number(atoms: ase.Atoms):
    return np.unique(atoms.get_atomic_numbers())

def get_nuclear_spin(num_n: int,num_p: int):
    return np.abs((num_n-num_p)/2)

def is_spin_active(nuclear_spin: float):
    return nuclear_spin!=0

def get_distance_from_point(positions: np.ndarray,ref: np.array=[0,0,0]):
    rel_positions=positions-ref
    return np.linalg.norm(rel_positions, axis=1)

#r_cutoff in angstrom
def within_r_cutoff(distance: np.array,r_cutoff: float=100):
    return distance <r_cutoff

