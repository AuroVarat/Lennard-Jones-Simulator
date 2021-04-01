import numpy as np

"""
CMod Ex1: PBC and MIC module.

Ascertains Periodic Boundary conditions and Minimum image cnvention on a simulation box and vector.

Author: Auro Varat Patnaik
Version: 11/03/2021
"""
def pbc(vector,box_size):
    """
        Constrains periodic boundary conditions

        :param vector: Any kind of vectorized element in list or on its own
        :param box_size: size of the simulation box

        :return PBC enforced vector: vector
        
    """
    return np.mod(vector,box_size)
def mic(vector,box_size):
    """
        Finds the image of the points with MIC

        :param vector: Any kind of vectorized element in list or on its own
        :param box_size: size of the simulation box

        :return MIC enforced vector image: vector
        
    """
    
    return np.subtract(np.mod(np.add(vector,np.multiply(box_size,0.5)),box_size),np.multiply(box_size,0.5))


