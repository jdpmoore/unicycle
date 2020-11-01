# AUTHOR
Sylvain Barbot (sbarbot@ntu.edu.sg)

# DESCRIPTION
Unicycle stands for **Unified Cycle of Earthquakes**. It computes fault slip and strain evolution
due to distant loading and coseismic stress change to model geodetic time series. Slip 
evolution is controlled by various spins of rate-and-state friction with radiation damping. 
Viscoelastic deformation is included using the analytic solution of

   Barbot S., J. D. P. Moore and V. Lambert, Displacement and Stress
   Associated with Distributed Anelastic Deformation in a Half Space,
   Bull. Seism. Soc. Am., 107(2), 10.1785/0120160237, 2017.

The method is efficient and accurate to simulate afterslip, slow-slip events, viscoelastic
relaxation and other aseismic processes.

# APPLICATIONS
* Forward calculation of fault slip evolution with
    * rate-strengthening friction, for afterslip calculation, including the effect of pre stress;
    * rate-and-state friction with radiation damping (slow-slip events).
* Stress and displacement kernels using the analytic solutions of
    * Okada (1985) and Okada (1992) for rectangular patches;
    * Meade (2007), Gimbutas et al. (2012), Nikkhoo and Walter (2015), for triangular patches.
    * Barbot, Moore, and Lambert (2017), for strain volumes.
* Kernel calculation using Andrew Bradley's hmmvp 1.3 for large matrices;
* Inversion of frictional parameters based on geodetic time series using simulated annealing;
* Prediction of remote stress, for calculation of Coulomb stress and other stress components;
* Inversion of residuals between geodetic time series and forward models;
* Include a list of 150 or so slip distributions of earthquakes and slow-slip events;
* Export input fault model (source and receivers) and simulation to Paraview in the legacy format.

# SOURCE
* Matlab implementation (OpenMP parallelism)
* Fortran90 implementation (MPI parallelism)

# TO DO
* simulation of InSAR data;
* accelerate computation with CUDAarray;
* neighborhood algorithm optimization;
* implement Linker and Dieterich rheology;
* incorporate pressure changes due to diffusion of surface precipitation (1D problem);
* add example plots of phase diagram with friction as a function of velocity;
* implement the state evolution of Kato and Tullis (2001);
* implement the ode23, ode113 solvers.
