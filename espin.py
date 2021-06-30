#!/usr/bin/env python
# coding: utf-8

# In[28]:


import numpy as np
from qutip import basis, spin_Jx, spin_Jy, spin_Jz


class espin:
    """
    A class to represent a central electron spin at a location.

    Attributes
    ----------
        mag_field:
            global magnetic field
        positions: 
            position of this spin. The coordinate should be on the same grid as that in the bath 
        gyromagnetic_ratio:
            gyromagnetic_ratio of the e spin. For now this is read in as we have experimental data
        state:
            spin state of the qubit. Currently initial state is set to the spin pointing in the x
            direction as a result of $\frac{\pi}{2}$ pulse
        hamiltonian:
            pure zeeman hamiltonian of the central spin. The assumption is that there is no 
            onsite spin-nuclear spin interaction because the spin is most likely on Carbon. 
    
    Usage Guide
        sp3=espin(B=np.array([0,0,5]))
        sp3.hamiltonian
    """
    def __init__(self, B, gyro: float = 2.0, position: np.array = [0,0,0]):
        """
        Constructs all the necessary attributes for the espin object.

        Parameters
        ----------
        mag_field:
            global magnetic field
        positions: 
            position of this spin. The coordinate should be on the same grid as that in the bath 
        gyromagnetic_ratio:
            gyromagnetic_ratio of the e spin. For now this is read in as we have experimental data
        """
        self.mag_field = np.array(B)
        self.gyromagnetic_ratio = gyro
        self.position = position
        self.state = (basis(2,0)+basis(2,1)).unit()
        self.hamiltonian = make_hamiltonian(self.gyromagnetic_ratio,self.mag_field)
        
def make_hamiltonian(gyro: float,B: np.array):
    """
        create pure zeeman hamiltonian. The assumption is that there is no 
        onsite spin-nuclear spin interaction because the spin is most likely on Carbon. 

        Parameters
        ----------
        B:
            global magnetic field
        gyro
            gyromagnetic_ratio of the e spin. For now this is read in as we have experimental data
        """
    S=np.array([spin_Jx(1/2), spin_Jy(1/2), spin_Jz(1/2)], 
        dtype=object)
    H_ms = gyro* B@S
    return H_ms

