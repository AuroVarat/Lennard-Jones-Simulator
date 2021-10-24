# Lennard-Jones Simulator

****All Parameters are in reduced units.**

## **A. File Structure:**

**Make sure to run the code from the directory containing all the three folders below.**

**1. Python:**

 - main.py: the main python script which is to be run
  
  -	Particle3D.py:library to create particles in 3D and manipulate them

 - Lennard.py : calculates force and potential under lennard jones potential using provided separation between particles

 - Mdutilities.py: generates a simulation box and assigns velocity and position to particles

-	PBC.py: contains methods to ensure periodic boundary conditions and minimum image convention

 **2. input-data:**

 - The folder contains. Txt files containing information about particles and external conditions. The naming convention of the files are <Type/Element>-<State>.txt e.g argon-solid.txt. The format of the text content discussed below.

 - argon-solid.txt and argon-gas.txt:contains initialization information for creating those particles.

 **3. output-data:**  
 
**The folder must be preset as the output data files are saved here.** 
Four types of files are saved here.  

 - trajectory.xyz : contains the trajectory of the particles for VMD
  visualization.

 - energy.dat : contains the time, potential energy, kinetic energy and
  total energy.

 - msd.dat : contains time and mean squared displacement as a function
  of time.
 - rdf.dat : contains radius and radial distribution average as a
  function of radius.

## **B. Quick Use:**

On calling the python file “main.py” , the source of input data must be provided, argon-solid and argon-gas are provided and saved in input-data folder.

Step 1:  
In terminal,

python main.py argon-solid

or if you have created your own particle-state follow the instructions in the next section for the format.

Step 2:

The main.py will run and save the data in the “output-data” folder.

## **C. How to create input files**

The output files must be .txt file.( This has been set for simplicity during command line argument and saving file name)

The input data file must contain:  
<Number of Particles> <Density*> <Tempreature> <time-step> <Number of Steps> <Cutoff>

<Name of Particle> <Mass of Particle*> **<x-position><y-position><z-position>< x-velocity><y-velocity ><z-velocity >** **The position and velocities can be random or just 0 as it is set in main.py using mdutilities.py.**
