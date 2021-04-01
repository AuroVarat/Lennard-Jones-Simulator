"""
CMod LJ: Simulation of Particles under Lennard Jones
Author: Auro Varat Patnaik
Version: 11/03/2021
"""
from Particle3D import Particle3D as p3d
import mdutilities 
import numpy as np
import sys
from lennard import lennard
#import cProfile



 
# Begin main code
def main():
    """ Read particle and paramters input file from command line"""
    
    if len(sys.argv)!=2:
        print("Wrong number of arguments.")
        print("Usage: " + sys.argv[0] + " <Input file Name>")
        quit()
    else:
        file_name = str(sys.argv[1])
        #itemp =str(sys.argv[2])

    """Set up simulation parameters"""
    p_handle = open("input-data/"+file_name+".txt","r") #reads input file 
    parameter = p_handle.readline().split() #reads first line containing parameters
    n = int(parameter[0]) #number of particle
    rho=float(parameter[1]) #Density
    temp=float(parameter[2]) #Tempreature
    dt = float(parameter[3]) #Time Step Size
    numstep=int(parameter[4]) # Number of Steps
    cutoff= float(parameter[5]) # cutoff radius
    time = 0.0 # initial time
    dr= 0.1 #for rdf dr
    #temp = itemp
    print("numstep ="+str(numstep))
    """Create Particles and set their initial conditions"""
    particle_list =p3d.new_particle(p_handle,n) #creates n particles 
    p_handle.close() #close input file
    box_size,_ = mdutilities.set_initial_positions(rho,particle_list) #set initial position of the particles and return box size
    mdutilities.set_initial_velocities(temp,particle_list) #set initial velocities
    separation_matrix,separation_modulus = p3d.get_separations(particle_list,box_size) #Find initial separation,force and potential energy
    force,p_e = lennard.get_interactions(separation_matrix,separation_modulus,cutoff) #Initial Forces and potential energy
    k_e = p3d.sys_kinetic(particle_list) #Find Initial Kinetic energy
    total_energy = k_e + p_e #Find initial total energy of the system 
    separation_modulus_alltime = np.zeros((numstep+1,n*n)) #Initialising an Array to contain all pair separation observed 
    separation_modulus_alltime[0] = separation_modulus #Adding initial pair separations
    #Saving the positions of initial state of system as reference for calculating MSD
    ref_state =[p.pos for p in particle_list]
    msd = p3d.get_msd(particle_list,ref_state,box_size) # calculating the initial MSD , expected to be zero
    
    """Initialise list to contain spatial( e.g RDF) and temporal data(e.g Energies over time)"""
    
    time_list = np.zeros(numstep+1)
    total_energy_list = np.zeros(numstep+1)
    potential_e_list = np.zeros(numstep+1)
    kinetic_e_list = np.zeros(numstep+1)
    msd_list = np.zeros(numstep+1)
    #Outfile file to record trajectory in VMD-compatible format
    outfile_pos = open("output-data/"+file_name+"-trajectory-"+str(dt)+".xyz", "w+")

    """Saving initial conditions"""
    time_list[0] = time
    total_energy_list[0] = total_energy
    potential_e_list[0] = p_e
    kinetic_e_list[0] = k_e
    msd_list[0] = msd
    outfile_pos.write("{0}\nPoint=1\n".format(n))  #Num Step of the system
    for particle in particle_list:
        outfile_pos.write(str(particle)+"\n")   #Position of all the particle in the nump step
 
    

    
    # Start the time integration loop
    for i in range(numstep):
        
        """ Velocity Verlet Integrator"""
        p3d.update_pos_2nd_list(particle_list,dt, force,box_size) # Update particle positions
        new_separation_matrix,new_separation_modulus = p3d.get_separations(particle_list,box_size) #calculate new separations
        force_new,p_e = lennard.get_interactions(new_separation_matrix,new_separation_modulus,cutoff) #calculate new forces
        p3d.update_vel_list(particle_list,dt, 0.5*(force+force_new)) #update velocities
        force = force_new # Re-define force value
        time += dt #Increase time
        "Setting and Saving temporal data to corresponding lists"
        #setting
        k_e = p3d.sys_kinetic(particle_list)
        msd = p3d.get_msd(particle_list,ref_state,box_size) 
        total_energy = k_e+p_e 
        
        
        #saving
        separation_modulus_alltime[i+1]= new_separation_modulus
        time_list[i+1] = time
        total_energy_list[i+1] = total_energy
        potential_e_list[i+1] = p_e
        kinetic_e_list[i+1] = k_e
        msd_list[i+1] = msd
        #For Trajectory
        outfile_pos.write("{0}\nPoint={1}\n".format(n,i+2))
        for particle in particle_list:
            outfile_pos.write(str(particle)+"\n")
      
    
    #Calculates the Radial Dsitrbution Average 
    g_r,radius_list = p3d.get_rdf(separation_modulus_alltime,box_size,dr,n,numstep,rho)
    #save spatial data(RDF) to file
    np.savetxt("output-data/"+file_name+"-rdf-"+str(dt)+".csv",np.c_[radius_list,g_r],fmt='%12.8f',header="radius,rdf")
    #Save Temporal data Lists to files
    np.savetxt("output-data/"+file_name+"-msd-"+str(dt)+".csv",np.c_[time_list,msd_list],fmt='%12.8f',header="time,msd")
    np.savetxt("output-data/"+file_name+"-energy-"+str(dt)+".csv",np.c_[time_list,potential_e_list,kinetic_e_list,total_energy_list],fmt='%12.8f',header="time,potential energy,kintic energy,total energy")
    
    # Close trajectory output file
    outfile_pos.close()
  
    
# Execute main method, but only when directly invoked
if __name__ == "__main__":
    main()
   #cProfile.run("main()")

