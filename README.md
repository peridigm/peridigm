# Peridigm [![CircleCI](https://circleci.com/gh/peridigm/peridigm.svg?style=shield)](https://circleci.com/gh/nperidigm/peridigm) [![Codacy Badge](https://app.codacy.com/project/badge/Grade/2d390b95b314447bac5673a07e7b4ea0)](https://www.codacy.com/gh/peridigm/peridigm/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=peridigm/peridigm&amp;utm_campaign=Badge_Grade)
[![Github Actions](https://github.com/peridigm/peridigm/actions/workflows/test.yml/badge.svg)](https://github.com/peridigm/peridigm/actions/workflows/test.yml)

Peridigm is an open-source computational peridynamics code developed, originally at Sandia National Laboratories and open-sourced in 2011, for massively-parallel multi-physics simulations.  It has been applied primarily to problems in solid mechanics involving pervasive material failure.  Peridigm is a C++ code utilizing foundational software components from Sandia's Trilinos project and is fully compatible with the Cubit mesh generator and Paraview visualization code.

An overview of Peridigm is given in [The Peridigm Meshfree Peridynamics Code](https://doi.org/10.1007/s42102-023-00100-0), by Littlewood, Parks, Foster, Mitchell, and Diehl, J Peridyn Nonlocal Model, 2023, and in the 2012 [Peridigm Users' Guide](http://www.sandia.gov/~djlittl/docs/PeridigmV1.0.0.pdf). Further details on software for computational peridynamics can be found in the [Roadmap for Peridynamic Software Implementation](http://www.sandia.gov/~djlittl/docs/PeridynamicSoftwareRoadmap.pdf).

Peridigm development began under the Physics & Engineering Models element of the US DOE's Advanced Simulation and Computing (ASC) program.  The project was led by Michael Parks and managed by John Aidun.  Subsequent funding has been provided by the US DOE through the ASC, ASCR, and LDRD programs.

## Examples
<table>
  <tr><td>
    <img width="1500" alt="Impact and Brittle Fracture" src="https://github.com/peridigm/peridigm/blob/master/doc/slide-image-1.jpg"/>
  </td><td>
    The simulation of impact and brittle fracture displayed here was achieved using explicit transient dynamics, the linear peridynamic solid constitutive model, short-range force contact, and a critical stretch bond failure law.  Peridynamics provides a natural framework for capturing pervasive material failure and fracture.
  </td></tr>
  <tr><td>
    <img width="1500" alt="Tensile Test Simulation" src="https://github.com/peridigm/peridigm/blob/master/doc/slide-image-2.jpg"/>
  </td><td>
    Peridigm is capable of performing explicit dynamic, implicit dynamic, and quasi-static time integration.  The tensile test simulation presented here was attained using an elastic correspondence constitutive model and quasi-static time integration.  Pre- and post-processing were carried out using Sandia's Cubit mesh generator and ParaView visualization code.
  </td></tr>
  <tr><td>
    <img width="1500" alt="Fragmentation of an Expanding Cylinder" src="https://github.com/peridigm/peridigm/blob/master/doc/slide-image-3.jpg"/>
  </td><td>
    The fragmentation of an expanding cylinder, shown here, was simulated using the linear peridynamic solid constitutive model and critical-stretch bond failure rule.  Initial velocities for each node in the discretization were specified via user-supplied analytic expressions.  Peridigm utilizes the RTCompiler function parser to process C-style expressions for the specification of input parameters, including intial and boundary conditions.
  </td></tr>
</table>

## Getting Started

Peridigm is a C++ code intended for use on Mac and Linux operating systems.  Both Peridigm and the Trilinos libraries it depends on should be built using MPI compilers and the [CMake](http://www.cmake.org/) build system.  The Peridigm test harness requires [python](https://www.python.org).  The build process has been tested using gcc and Intel compilers, Open MPI, and MPICH.  The steps below should be followed in order, beginning with installation of the required third-party libraries.

 * [Installing Third-Party Packages and Libraries](https://github.com/peridigm/peridigm/blob/master/doc/InstallingThirdPartyLibs.md)
 * [Building Peridigm](https://github.com/peridigm/peridigm/blob/master/doc/BuildingPeridigm.md)
 * [Running Simulations with Peridigm](https://github.com/peridigm/peridigm/blob/master/doc/RunningSimulations.md)
 * [Running Simulations with Peridigm Docker Image](https://github.com/peridigm/peridigm/blob/master/doc/RunningSimulationsDocker.md)


## Team

The following individuals have made significant contributions to the Peridigm code:
*  [Michael Parks](https://cfwebprod.sandia.gov/cfdocs/CompResearch/templates/insert/profile.cfm?snl_id=109440&ename=Parks,+Michael+L.) ([@mlparks](https://github.com/mlparks))
*  [David Littlewood](https://cfwebprod.sandia.gov/cfdocs/CompResearch/templates/insert/profile.cfm?snl_id=195431&ename=Littlewood,+David+John) ([@djlittl](https://github.com/jdlittl))
*  [John Mitchell](https://cfwebprod.sandia.gov/cfdocs/CompResearch/templates/insert/profile.cfm?snl_id=13850&ename=Mitchell,+John+A.) ([@jamitch-snl](https://github.com/jamitch-snl))
*  [Stewart Silling](http://www.sandia.gov/~sasilli/)
*  [John Aidun](http://www.cs.sandia.gov/materials_methods/)
*  [Daniel Turner](https://cfwebprod.sandia.gov/cfdocs/CompResearch/templates/insert/profile.cfm?snl_id=121696&ename=Turner,+Daniel+Z.)
*  [John Foster](http://johnfoster.pge.utexas.edu/) ([@johntfoster](https://github.com/johntfoster))

## Citing

When citing Peridigm, please reference the following:

M.L. Parks, D.J. Littlewood, J.A. Mitchell, and S.A. Silling, Peridigm Users' Guide, Tech. Report SAND2012-7800, Sandia National Laboratories, 2012.

## License

See [our license information](https://github.com/peridigm/peridigm/blob/master/LICENSE.md).
