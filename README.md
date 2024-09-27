# IGAMDS
Isogeometric solver for the monodomain model (electrophysiological reaction-diffusion system) over NURBS surfaces. The monodomain equation is integrated in space via a NURBS-based isogeometric method, in time via the semi-implicit Crack-Nicholson/Admas-Bashforth scheme. Derivative over the surface are computed using the local curvilinear grame defined by the metric tensor. The associated ionic current model is advanced in time via a fourth-order Runge-Kutta scheme. Three possible ionic current models are available: Aliev-Panfilov, Rogers-McCulloch, Beeler-Reuter.

# Developers
Alessandro Nitti, Polytechnic University of Bari (https://scholar.google.it/citations?user=lv1V6-4AAAAJ&hl=it&oi=ao)  

# Methodology
The governing equations are the following: 

$$ C_m \frac{\partial v}{\partial t} - \nabla \cdot \left( D \nabla v \right) + \chi i_{ion} = \chi i_a \, $$ 

Details about the methodology and the implementation can be found in the following papers:  
https://doi.org/10.1016/j.cma.2021.113877 
https://doi.org/10.1016/j.cma.2004.10.008 

# How to run a test
1. Create the desired NURBS surface and performs both h-refinement (element insertion) and p-refinement (degree elevation of basis functions) by running the script ./proproc/MAIN_<test_name>.m. This will generate the test folder with the related input files.
2. Insert the desired values of physical parameters and select the desired ionic current model via <test_name>/input/modelpar.in.
3. Prescribe the stimulation protocol by means of the file <test_name>/input/stim_prot.in. 
4. Compile the ./src/ files via makefile  
5. Run the test in the <test_name> folder

# Organization of the repository
./preproc/: contains matlab pre-processing scripts that generate the test folder and the input files needed for the execution. It is used to design the Eulerian and Lagrangian grid.  
./postproc/: contains matlab post-processing scripts that allow to visualize the simulation output. Pressure/velocity/vorticity fields are displayed along with the body position.  
./src/: contains the makefile and all the routines to execute the program.  
./test_planeslab/: contains the input files for the test case simulating the propagation of the action potential over a plane slab with the Aliev-Panfilov ionic current model. The time-traces of the action potential is verified against data provided in referece https://doi.org/10.1007/s00466-009-0434-z.
./test_spiral/: contains the input files for the test case simulating the trigger of a self-sustained spiral action potential. The Aliev-Panfilov ionic current model is used.
./test_cyl/: contains the input files for the test case simulating the propagation of an Aliev-Panfilov action potential over a cylindrical shell. The anisotropic diffusion occurs as a result of the surface geometry.
