from qutip import basis
class central_spin(positions: np.array, ms: float):
   self.positions = positions
   self.ms = ms
   #initialize along x
   self.state_vector = (basis(2,0)+basis(2,1)).unit()

