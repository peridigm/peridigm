# Peridigm

Peridigm is an open-source computational peridynamics code developed at Sandia National Laboratories for massively-parallel multi-physics simulations.  It has been applied primarily to problems in solid mechanics involving pervasive material failure.  Peridigm is a C++ code utilizing foundational software components from Sandia's Trilinos project and is fully compatible with the Cubit mesh generator and Paraview visualization code.

The [peridigm-users](https://software.sandia.gov/mailman/listinfo/peridigm-users) e-mail list connects Peridigm enthusiasts and provides a forum for user questions.

The 2012 [Peridigm Users' Guide](http://www.sandia.gov/~djlittl/docs/PeridigmV1.0.0.pdf) gives an overview of Peridigm's core capabilities. Further details on software for computational peridynamics can be found in the [Roadmap for Peridynamic Software Implementation](http://www.sandia.gov/~djlittl/docs/PeridynamicSoftwareRoadmap.pdf).

Peridigm development began under the Physics & Engineering Models element of the US DOE's Advanced Simulation and Computing (ASC) program.  The project was led by Michael Parks and managed by John Aidun.  Subsequent funding has been provided by the US DOE through the ASC, ASCR, and LDRD programs.


## Getting Started

Peridigm is a C++ code intended for use on Mac and Linux operating systems.  Both Peridigm and the Trilinos libraries it depends on should be built using MPI compilers and the [CMake](http://www.cmake.org/) build system.  The Peridigm test harness requires python.  The build process has been tested using gcc and Intel compilers, Open MPI, and MPICH.  The steps below should be followed in order, beginning with installation of the required third-party libraries.

 * [Installing Third-Party Packages and Libraries](https://github.com/peridigm/peridigm/blob/master/doc/InstallingThirdPartyLibs.md)

 * [Building Peridigm](https://github.com/peridigm/peridigm/blob/master/doc/BuildingPeridigm.md)

 * [Running Simulations with Peridigm](https://github.com/peridigm/peridigm/blob/master/doc/RunningSimulations.md)

