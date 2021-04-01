"""

 CompMod Ex2: Particle3D, a class to describe point particles in 3D space

 An instance describes a particle in Euclidean 3D space: 
 velocity and position are [3] arrays

 Includes time integrator methods, dynamics methods for individual particle or list of particles.

Author: Auro Varat Patnaik
Version: 11/03/2021

"""

import numpy as np
import PBC
# Constants should go here
G0 = 6.67384E-11            # From CODATA 2010
ASTRO_U = 149597870700.0    # From IAU resolution 2012/08 B2
YEAR = float(3.15576e7)     # Julian year = 365.25 days


class Particle3D(object):
    """
    Class to describe point-particles in 3D space

        Properties:
    label: name of the particle
    mass: mass of the particle
    pos: position of the particle
    vel: velocity of the particle

        Methods:
    __init__
    __str__
    kinetic_e  - computes the kinetic energy
    momentum - computes the linear momentum
    update_pos_1st - updates the position to 1st order
    update_pos_2nd - updates the position to 2nd order
    update_vel - updates the velocity

        Static Methods:
    new_particle - initializes a P3D instance from a file handle
    sys_kinetic - computes total K.E. of a p3d list
    com_velocity - computes total mass and CoM velocity of a p3d list
    """

    def __init__(self, label, mass, pos, vel):
        """
        Initialises a particle in 3D space

        :param label: String w/ the name of the particle
        :param mass: float, mass of the particle
        :param position: [3] float array w/ position
        :param velocity: [3] float array w/ velocity
        """

        self.label = str(label)
        self.mass = float(mass)
        self.pos = np.array(pos,float)
        self.vel = np.array(vel,float)


    @staticmethod
    def new_particle(file_handle,n=1):
        """
        Initialises a Particle3D instance given an input file handle.
        
        The input file should contain one per particle in the following format:
        label   <mass>  <x> <y> <z>    <vx> <vy> <vz>
        
        :param input file path: str, path to readable file containting data in the above format
        :param n: int, Number of particles required
        
        :return: Particle3D instance/ List of Particle3D instance
        """
        
        #Get File Handle from path
        particle_info = file_handle.readline().split()
        #Takes data from the file line
        particle_label = particle_info[0]
        particle_mass = particle_info[1]
        particle_position = particle_info[2:5]
        particle_velocity = particle_info[5:8]
        
        if n == 1 :
            return Particle3D(particle_label,particle_mass,particle_position,particle_velocity)
        elif n > 1:
            list_particles = []
            for i in range(n):
                list_particles.append(Particle3D(particle_label+str(i),particle_mass,particle_position,particle_velocity))
            return list_particles
        elif n == 0 :
            return Particle3D("Tachyon",0.,np.zeros(3),np.zeros(3))  #creates a tachyon   


    def __str__(self):
        """
        XYZ-compliant string. The format is
        <label>    <x>  <y>  <z>
        """
         
        return str(self.label)+" "+str(self.pos[0])+" "+str(self.pos[1])+" "+str(self.pos[2])


    def kinetic_e(self):
        """
        Returns the kinetic energy of a Particle3D instance

        :return ke: float, 1/2 m v**2
        """
        ke = 0.5*self.mass*(np.linalg.norm(self.vel)**2)
        return ke


    def momentum(self):
        """
        Returns the linear momentum of a Particle3D instance
        :return p: (3) float np.array, m*v
        """
        p = self.mass*self.vel
        return p


    def update_pos_1st(self, dt):
        """
        1st order position update

        :param dt: timestep
        """
        self.pos = self.pos + self.vel*dt
        
      


    def update_pos_2nd(self, dt, force):
        """
        2nd order position update

        :param dt: timestep
        :param force: [3] float array, the total force acting on the particle
        """

        self.pos = self.pos + self.vel*dt + 0.5*(1/self.mass)*force*(dt**2)
        
    @staticmethod
    def update_pos_2nd_list(p3d_list, dt, force_list,box_size):
        """
        2nd order position update

        :param dt: timestep
        :param force: [3] float array, the total force acting on the particle
    
        """
       
        for p3d,force in zip(p3d_list,force_list):
            p3d.pos = PBC.pbc(p3d.pos +p3d.vel*dt + 0.5*(1/p3d.mass)*force*(dt**2) ,box_size)
       
    

    def update_vel(self,dt, force):
        """
        Velocity update

        :param dt: timestep
        :param force: [3] float array, the total force acting on the particle
        """
     
        self.vel = self.vel + (1/self.mass)* force*dt
        
    @staticmethod
    def update_vel_list(p3d_list,dt, force_list):
        """
        Velocity update

        :param dt: timestep
        :param force: [3] float array, the total force acting on the particle
        """
        n = len(p3d_list)
        
        for p3d,force in zip(p3d_list,force_list):
            
            p3d.vel = (p3d.vel + (1/p3d.mass)* force*dt)

    @staticmethod
    def sys_kinetic(p3d_list):
        """
        Returns the kinetic energy of the whole system

        :param p3d_list: list in which each item is a P3D instance

        :return sys_ke: \sum 1/2 m_i v_i^2 
        """
        sys_ke = 0
        for p3d in p3d_list:
            sys_ke += p3d.kinetic_e()
        return sys_ke


    @staticmethod
    def com_velocity(p3d_list, with_mass = False):
        """
        Computes the total mass and CoM velocity of a list of P3D's

        :param p3d_list: list in which each item is a P3D instance
        :param with_mass: boolean, optional. 
           Switch determining nature of return value. When it is False (default) just the Centre-of-mass velocity
           is returned, when True total mass of the particles in the list is also returned.

        :return com_vel: float, Centre-of-mass velocity
        :return total_mass: float, The total mass of the system
            Present in a tuple with com_vel only if with_mass = True.
        """
        total_mass = 0
        total_momentum = 0


        for p3d in p3d_list:
            total_mass += p3d.mass
            total_momentum += p3d.momentum()
        com_vel = total_momentum/total_mass

        if with_mass == True:
            return com_vel,total_mass
        else:
            return com_vel
        
    @staticmethod
    def get_msd(p3d_list,ref_list,box_size):
        """
        Computes the Mean Squared Displacement of Particles

        :param p3d_list: list in which each item is a P3D instance
        :param ref_list: list in which each item is a Position of the initial state of P3D instances

        :return MSD: float, MSD of the time
        
        """
        dis=0
        for i,o in zip(p3d_list,ref_list):
            dis += np.square( np.linalg.norm(  PBC.mic(i.pos-o,box_size)   )   )
        return dis/len(p3d_list)
    
    
    
    @staticmethod
    def get_rdf(separations_list,box_size,dr_element,number_of_particles, numstep,rho):
        """
        Computes the Radial Distribution Average of Particles

        :param separations_list: list in which each item is a modulus of separation of P3D instances
        :param box_size: Size of the simulation cube
        :param dr_element: Size of the small radial element for Density function
        :param number_of_particles: Total number of particles in the simulation
        :param numstep: Total runs
        :param rho: Density of the system
        

        :return g(r): float array,  Normalized frequency of occurence of separations
        :return r: float array,  list of occuring Radius of separations 
        
        """
        bins = np.arange(dr_element,(box_size[0]*(3**(0.5))/2)+dr_element,dr_element) #creatings bins from minimum radial length to sphere contaning the simulation cube/diagonal of the cube
        hist,_ = np.histogram(separations_list[separations_list!=0],bins=bins,density=False)
        r_list = bins[:-1]
        g_r = Particle3D.normalize_hist(hist,r_list,dr_element,number_of_particles,numstep,rho)
        return g_r,r_list
    
    @staticmethod 
    def normalize_hist(histogram_data,radius_list,dr_element,number_of_particles, numstep,rho):
        """
        Normalizes histogram data / Time Averaged pair separations

        :param histogram_data: list of frequencies of occurence of pair separations
        :param radius_list: list of radius of separatoins in considertation
        :param box_size: Size of the simulation cube
        :param dr_element: Size of the small radial element for Density function
        :param number_of_particles: Total number of particles in the simulation
        :param numstep: Total runs
        :param rho: Density of the system
        

        :return g(r): float array,  Normalized frequency of occurence of separations
      
        
        """ 
        density = 4*np.pi*rho*(np.square(radius_list+dr_element*0.5))*dr_element #denisty of the ideal homogenous system
        g_r = (2*(histogram_data))/(density*number_of_particles*(numstep+1)) #Normalized histogram
        return g_r
    
    @staticmethod
    def get_separations(p3d_list,box_size):
        """
        Computes the Pair Separations present in between particles of the system

        :param p3d_list: list in which each item is a P3D instance
        :param box_size: float array, size of the simulation box

        :return separation matrix: N x N x 3 array, containing pair separations
        :return separation modulus matrix: N x N  array, containing modulus of pair separations
        
        """
        N=len(p3d_list)
        p3d_pos_list = np.zeros((N,3))
        #Creates a position list
        for index,p3d in enumerate(p3d_list):
            p3d_pos_list[index] = p3d.pos
        separation_matrix = np.ones((N,N,3))*p3d_pos_list
        #Calculates all pair separations and ascetains MIC
        separation_matrix = PBC.mic(np.swapaxes((p3d_pos_list-np.swapaxes(separation_matrix, 0, 1)), 0, 1),box_size )
        #Removes the lower diagonal triangle of the matrix - This will enable us to reduce number of calculations by taking advantage of the symmetry
        separation_matrix[np.triu_indices(n=N, k=0)]=0
        separation_modulus_matrix = np.linalg.norm(separation_matrix,axis=2)
    
        return separation_matrix,separation_modulus_matrix.flatten()
  