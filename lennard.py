
"""
CMod Ex3: Force and Potential Model for Morse Potential.
Author: Auro Varat Patnaik
Version: 11/03/2021
"""

import numpy as np


class lennard():
    
    @staticmethod
    def interaction(separation_vec,separation_mod,cutoff):
        """
        Method to return the force and Potential on a particle
        in a Lennard Jones potential.
        Force is given by
        F(x) = -dV/dx = 48*(1/r**14 - 1/2*r**8)(r_vec)

        :param separation_vec: Pair Separation independant or Array
        :param separation_mod: Modulus of pair Separation independant or Array
        :param cutoff: Cutoff of the system 
        
       
        :return: force and potential acting on particle as Numpy array
        """
        separation_mod_for_force = np.where(separation_mod < [cutoff],separation_mod,0)
        #New Line
        separation_mod_for_potential = np.where(separation_mod < [cutoff],separation_mod,[cutoff])
        #turns all separations beyond cutoff to negligent
        #separating common terms for easier/quicker calculatoins
        
        s2 = np.divide(1,np.square(separation_mod_for_force),
                       out = np.zeros_like(separation_mod_for_force),where=separation_mod_for_force!=0)
        s6 = np.power(s2,3)
        #New Line
        s6_for_pot = np.divide(1,np.power(separation_mod_for_potential,6),
                       out = np.zeros_like(separation_mod_for_potential),
                       where=separation_mod_for_potential!=0) 
      
        potential = 4*(np.square(s6_for_pot)-s6_for_pot)
        force = 48*(s2*np.square(s6)-0.5*s6*s2)*separation_vec 
        
        
        return force,potential 

    @staticmethod
    def get_interactions(separation_matrix,separation_modulus_matrix,cutoff):
        """
        Method to return the force and potential on a list of particle under Lennard Jones Potential
        

        :param separation_matrix: N x N x 3 matrix containing all Pair Separations between particles
        :param separation_modulus: N x N matrix containing all Modulus of pair Separations 
        :param cutoff: Cutoff of the system 
        
       
        :return: force and potential acting on all particle as Numpy array
        """
        N = len(separation_matrix)
        force_Matrix,potential_matrix = lennard.interaction(separation_matrix.reshape(N*N,3),separation_modulus_matrix.reshape(N*N,1),cutoff)    
        force_Matrix = force_Matrix.reshape(N,N,3) #Reshapes to introduce particular particle informatoin
        force_Matrix -= np.swapaxes(force_Matrix, 0, 1)  
        force_list = np.sum(force_Matrix,axis=1)
        potential = np.sum(potential_matrix) #sums to the total potential energy of the system
        
        return force_list,potential
            
        
    