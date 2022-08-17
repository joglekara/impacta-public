# impacta
```impacta``` (Implicit Magnetized Plasma And Collisional Transport with Anisotropy) is a fully-implicit code for electron heat transport in inertial fusion relative scenarios. It solves the Vlasov-Fokker-Planck equations under the Lorentz approximation for electrons using a Cartesian Tensor expansion up to second order, with third order added in a reduced model. Additional physics models include a simple hydrodynamic solver for ions, ionization and a ray tracing package for laser heating using inverse bremsstrahlung.

## Contributors:
The code was based on ```impact```, conceived and developed by Robert Kingham and Tony Bell (R.J. Kingham and A.R. Bell, J. Comp. Phys. 2004).
A new C++ version of the code with the addition of second order Cartesian Tensor terms was developed by Alec Thomas.
Hydrodynamic ion motion were added by Alec Thomas based on a model developed for ```impact``` by Chris Ridgers.
A ray tracing package was added by Archis Joglekar.

## Installation instructions:
The code requires PETSc and boost libraries
