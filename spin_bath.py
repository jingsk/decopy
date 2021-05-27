from bath.py import bath
import numpy as np
from qutip import basis

class spin_bath(bath):
   """
   makes bath hamiltonian based on cce order
   should be created after knowing cce order 
   (also for converging cce_order)
   """
   def __init__(self, gyromagneticratio_data_set)
   self.bath_sites=bath
   self.gyro=gyro_data_set
   self.state_vector=generate_initial_state(bath_site.get_number_of_sites())
   self.hamiltonian=np.array([[np.eye(2)]*bath_site.get_number_of_sites()])

   @classmethod
   def create_from_cce(cls, cce_order):
      self.hamiltonian=make_hamiltonian(bath_sites.get_positions(), cce_order)

def generate_initial_state(n_sites):
   """
   randomly assign up or down magnetic state
   """
   state_list=[]
   for site in range(n_site):
      state=basis(2, intnp.rand())
      state_list.append(state)
   return np.array(state_list)

def make_hamiltonian(positions, cce_order):
   """
   for cce=2, for each site, loop thru other sites and make h_nn
   """
   h_list=[]
   if cce_order ==2:
      h_list=[]
      for site in range(positions.shape/3):
         for other_site in range(positions.shape/3):
               h_list.append(np.eye(2)+get_nn_interactions(site1, site2))
      return np.array(h_list)

def get_nn_interactions(site1, site2, gyro):
   """
   calculate h_nn based on eq in paper 
   """
   r1=site1.positions
   r2=site2.positions
   I1=state_vector(site1.index)
   T2=state_vector(site2.index)
   #add equation to calculate interaction hamiltonian
   return site_hamiltonian

