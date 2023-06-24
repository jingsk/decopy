from bath.py import bath
import numpy as np
from qutip import basis

class spin_bath(bath):
    """
    makes bath hamiltonian of each bath site.
    make_hamiltonian should be used after knowing 
    cce order perhaps called by CCE class.
    """
    def __init__(self, bath, gyromagneticratio_data_set):
        self.bath_sites=bath
        self.gyro=gyro_data_set
        self.state_vector=generate_initial_state(bath_site.get_number_of_sites())
        self.hamiltonian=np.array([[np.eye(2)]*self.bath_site.get_number_of_sites()])

def generate_initial_state(n_sites):
    """
    randomly assign up or down magnetic state
    """
    state_list=[]
    for site in range(n_site):
        state=basis(2, int(np.rand()))
        state_list.append(state)
    return np.array(state_list)

def make_hamiltonian(positions, cce_order, bath_vertices):
    """
    for cce=2, for each site, loop thru other sites and make h_nn
    """
    h_list=[]
    vertices = bath_vertices[cce_order]
    if cce_order ==2:
        h_list=[]
        for site in np.unique(vertices[:,0]):
            array_with_site = a[np.where(vertices[:,0]==site)]
            for other_site in array_with_site[0,:]:
                h_list.append(np.eye(2)+get_nn_interactions(site, other_site))
        return np.array(h_list)

def get_nn_interactions(site1, site2, gyro):
    """
    calculate h_nn based on eq in paper 
    currently not very complete
    """
    r1=site1.positions
    r2=site2.positions
    I1=state_vector(site1.index)
    I2=state_vector(site2.index)
    #add equation to calculate interaction hamiltonian
    return site_hamiltonian
