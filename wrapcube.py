import numpy as np
from ase.io.cube import read_cube_data, write_cube
from ase.geometry import wrap_positions
from ase.visualize import view

#the wrap part was modified from wrap_positions in ase.geometry
def wrap_and_center(positions, cell, pbc=np.array([True,True,True]), nxnynz=[100,100,100],eps=1e-7):
    shift = np.zeros(3) - eps
    # Don't change coordinates when pbc is False
    shift[np.logical_not(pbc)] = 0.0
    #positions in crystal coordinate
    fractional = np.linalg.solve(cell.T,np.asarray(positions).T).T - shift
    #for translation later
    nxnynz_moved=np.zeros(3,dtype=int)
    #for each crystal axis
    for i in range(3):
        if not pbc[i]:
            continue
        # translation increment for both atoms and volume
        increment=1/nxnynz[i]
        
        #wrap part
        #indices for reordering below
        indices = np.argsort(fractional[:, i])
        #reordered fragment x,y,z position
        sp = fractional[indices, i]
        #take the difference between the largest and the smallest
        widths = (np.roll(sp, 1) - sp) % 1.0
        #minus the smallest of the differences
        diff=sp[np.argmin(widths)]
        #total translation once corrected for increment
        diff_inc=(np.floor_divide(diff,increment))*increment
        nxnynz_moved[i]=diff_inc/increment
        #move
        fractional[:, i] -= diff_inc
        #wrap
        fractional[:, i] %= 1.0

        #now move to center
        current_center=(np.max(fractional[:, i])-np.min(fractional[:, i]))/2
        desired_center=0.5
        translate=desired_center-current_center
        translate_inc=(np.floor_divide(translate,increment))*increment
        fractional[:, i] += translate_inc
        nxnynz_moved[i]-=translate_inc/increment
        nxnynz_moved[i]=nxnynz_moved[i] % nxnynz[i]
    return np.dot(fractional, cell), nxnynz_moved

#for reading in geometry and density from a cube file, wraping, centering, and writing 
def read_wrap_write(path, file):
    density, atoms = read_cube_data(path+file)
    nxnynz=np.shape(density)
    pos=atoms.get_positions()
    translate_directions=np.array([True,True,False])
    pos2,nxnynz_moved=wrap_and_center(atoms.get_positions(),atoms.get_cell(),pbc=translate_directions,nxnynz=nxnynz)
    atoms.set_positions(pos2)
    #view(atoms) for debugging
    #roll density forward with a negative number
    density2=np.roll(density,-nxnynz_moved,axis=[0,1,2])
    file2=file.split('.')[0]+'_centered.cube'
    with open(path+file2,'w') as f:
        write_cube(f,atoms,data=density2)

#a useful function for exporting 
def cleanUp(atoms):
    new_position, dummy = wrap_and_center(atoms.get_positions(),atoms.get_cell(),pbc=atoms.get_pbc())
    atoms.set_positions(new_position)
    return atoms
#bs='/projectnb/fpmats/jing/CNT/bandstructure/'
#pathFileList=[
    #(bs+'CNTligand/CNTligandOptimizeLDA_col_SCF_gammaonly/split/','CHGCAR_mag.cube')
    #(bs+'CNTligand1/CNTligand1LDAcol_SCF/split/','CHGCAR_mag.cube'),
    #(bs+'CNTligand2/CNTligand2LDAcol_SCF/split/', 'CHGCAR_mag.cube'),
    #(bs+'CNTligandBr/CNTligandBr_LDA_col_SCF/split/','CHGCAR_mag.cube'),
    #(bs+'CNTligandCOOH/CNTligandCOOH_LDA_col_SCF/split/', 'CHGCAR_mag.cube'),
    #(bs+'CNTligandNO2/CNTligandNO2_LDA_col_SCF/split/', 'CHGCAR_mag.cube')
    #]
#for path, file in pathFileList:
#    read_wrap_write(path,file)