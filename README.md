# decopy
calculate T2 in a system using the spin-diffusion model. We had started implementing our own CCE code which produces the coherence function given a structure. Our Python package integrated with the atomic simulation environment (ASE) to quickly return geometric parameters and natural isotropic abundance data. To describe the spin hamiltonian in the basis of electron, nuclei spins we use the  Quantum Toolbox in Python package (QuTiP) package (https://github.com/qutip/qutip) which can handle basis, spin operators generation. QuTiP runs on a Cython engine which runs at C-code speed making it an attractive choice for linear algebra calculations. In 2021, Mykyta Onizhuk and Giulia Galli released their complete code, PyCCE,\cite{Onizhuk2021a} which led us to stop development of our code. 

## file descriptions
espin.py - central spin class - define electron spin Hamiltonian (ZFS, magnetic field splitting) ($H_s$).
spinbath.py - bath spin class - define nuclear bath Hamiltonian (ZFS, magnetic field splitting, spin-spin interactions) ($H_b$).
bath.py - bath class - reads in geometry file (rn CONTCAR), randomly assign isotope, retain only spin active site from a ref point with a cutoff radius.
util:
wrapcube.py - wrap and center code originally written for wrapping density in .cube files. 

# Theoretical background
## Spin Hamiltonian of a central spin
To account for interactions among the central spin and all the spin-active nuclei spins, we write out the interactions in the following way:
To account for interactions among the central spin and all the spin-active nuclei spins, we write out the interactions in the following way:
$$\hat{H} = \hat{H}\_S +\hat{H}\_{SB}+\hat{H}\_B$$,    
where $\hat{H}\_S$, $\hat{H}\_{SB}$, and $\hat{H}\_B$ are the spin Hamiltonians associated with the spin, the spin-bath interactions, and the bath. They are defined as follows
$$\hat{H}\_S = \mathbf{SDS}+\mathbf{B}\gamma\_S \mathbf{S}$$,
$$\hat{H}\_{SB} = \sum\_{I}\mathbf{S}\mathbf{A}\_i\mathbf{I}\_i$$,
$$\hat{H}\_{B} = \sum\_{i}\mathbf{I}\_i\mathbf{P}\_i\mathbf{I}\_i +\mathbf{B}\gamma\_i \mathbf{I}\_i+\sum\_{I>j}\mathbf{I}\_i\mathbf{J}\_{ij}\mathbf{I}\_j$$,
$$\mathbf{J}\_{ij}=\frac{\mu\_0}{4\pi} \gamma\_i \gamma\_j \hbar^2\left ( \frac{\mathbf{\vec{I}}\_i\cdot \mathbf{\vec{I}}\_j}{|r\_{ij|^3}}-\frac{3(\mathbf{\vec{I}}\_i\cdot \vec{r}\_{ij})(\mathbf{\vec{I}}\_j\cdot \vec{r}\_{ij})}{|r\_{ij}|^5} \right )$$,
with $\mathbf{S} = (\hat{S}\_x,\hat{S}\_y,\hat{S}\_z)$ and $\mathbf{I} = (\hat{I}\_x,\hat{I}\_y,\hat{I}\_z)$,the components of the spin operators of the central and a bath nuclei spin, $\gamma\_S$ and $\gamma\_i$, the gyromagnetic ratio of the central spin and of the nuclei index i., $A\_i$, the hyperfine tensor between the central spin and the nuclei index I, $D$($P\_i$), the self-interaction terms of the central spin, quadrupole tensor of nuclei i, $J\_{ij}$, the interaction tensors between bath spins. Instead of calculating the hyperfine tensors explicitly in DFT for all nuclei, we use the point-dipole approximation to estimate the interactions. This assumption is valid for localized electron spins. Employing the point-dipole approximation,
$$\mathbf{A}\_{i}=\frac{\mu\_0}{4\pi} \gamma\_i \gamma\_s \hbar^2\left ( \frac{\mathbf{I}\_3}{|r\_{i}|^3}-\frac{3(\vec{r}\_{i}\otimes \vec{r}\_{i})}{|r\_{i}|^5} \right ),$$

where $\mathbf{I}$ is a 3x3 identity matrix.

It is possible to project the hamiltonian in the basis of the two-level electron spins simplifying the hamiltonian making $\mathbf{SDS}$ the zero-field splitting (ZPS), which can be obtained from experiment. 
