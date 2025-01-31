# IGAMDS
Isogeometric solver for the monodomain model (electrophysiological reaction-diffusion system) over NURBS surfaces. The monodomain equation is integrated in space via a NURBS-based isogeometric method, in time via the semi-implicit Crank-Nicholson/Admas-Bashforth scheme. Derivatives over the surface are computed using the local curvilinear frame defined by the metric tensor. The associated ionic current model is advanced in time via a fourth-order Runge-Kutta scheme. Three possible ionic current models are available: Aliev-Panfilov (Aliev, R. R., & Panfilov, A. V. (1996). A simple two-variable model of cardiac excitation. Chaos, Solitons & Fractals, 7(3), 293-301.), Rogers-McCulloch (Rogers, J. M., & McCulloch, A. D. (1994). A collocation-Galerkin finite element model of cardiac action potential propagation. IEEE Transactions on Biomedical Engineering, 41(8), 743-757.), Beeler-Reuter (Beeler, G. W., & Reuter, H. (1977). Reconstruction of the action potential of ventricular myocardial fibres. The Journal of physiology, 268(1), 177-210.).

# Developers
Alessandro Nitti, Polytechnic University of Bari (https://scholar.google.it/citations?user=lv1V6-4AAAAJ&hl=it&oi=ao)  

# Methodology
The governing equations are the following: 

$$ 
\begin{aligned}
& C_m \frac{\partial v}{\partial t} - \nabla \cdot \left( \mathbf{D} \nabla v \right) + \chi \ i_{ion} = \chi \ i_a \\
& i_{ion}=f(v,w) \\
& \frac{dw}{dt}=g(v,w) 
\end{aligned}
$$  

$$
\begin{aligned}
& x \rightarrow \text{[cm], space}  \\
& t \rightarrow \text{[ms], time} \\
& v \rightarrow \text{[mV], action potential} \\
& w \rightarrow \text{[-], recovery variable} \\
& C_m  \rightarrow \text{[mF/cm}^3 \text{], tissue capacity } \\
& \mathbf{D}  \rightarrow \text{[mA/(mV*cm)], conductivity tensor } \\
& \chi  \rightarrow \text{[1/cm], surface to volume ratio} \\
& i_{ion}  \rightarrow \text{[mA/cm}^2 \text{], ionic current density} \\
& i_a \rightarrow \text{ [mA/cm} ^2 \text{], applied current density}
\end{aligned}
$$ 

Details about the methodology and the implementation can be found in the following research papers:  
https://doi.org/10.1016/j.cma.2021.113877  
https://doi.org/10.1016/j.cma.2004.10.008  

The conducticity tensor is provided in the parametric space, i.e., on Cartesian bases, and than transformed to the local curvilinear bases on the surface.

# Dependencies 
Install the HDF5 library (sudo apt get ), any release from 2015 on.  
Linear systems are solved via the MGMRES algorithm provided by John Burkardt, (Barrett R. et al., Templates for the Solution of Linear Systems: Building Blocks for Iterative Methods, SIAM, 1994. ISBN: 0898714710, LC: QA297.8.T45.)
   
# How to run a test
1. Create the desired NURBS surface and performs both h-refinement (element insertion) and p-refinement (degree elevation of basis functions) by running the script ./proproc/MAIN_<test_name>.m. This will generate the test folder with the related input files.  
2. Prescibe the desired physical parameters, select the ionic current model and prescribe the stimulation protocol in the write-to-file section of the ./proproc/MAIN_<test_name>.m file.
5. Compile the ./src/ files via makefile    
6. Run the test in the <test_name> folder  

# Organization of the repository
./preproc/: contains matlab pre-processing scripts that generate the test folder and the input files needed for the execution. It is used to design the Eulerian and Lagrangian grid.    
./postproc/: contains matlab post-processing scripts that allow to visualize the simulation output. Pressure/velocity/vorticity fields are displayed along with the body position.   
./src/: contains the makefile and all the routines to execute the program.  
./test_cyl.mp4: tape showing the action potential propagation on a cylindrical shell.
./test_spiral.mp4: tape showing the trigger of a self-sustained spiral wave potential.
./test_planeslab.png: verification of the action potential time traces against literature data.
