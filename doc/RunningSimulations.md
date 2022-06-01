# Running Simulations with Peridigm

Peridigm is run from the command line and requires both an input deck (text file) and a discretization.  Input decks may be in either `.yaml` format or `.xml` format.  The content of both file formats is identical, the main difference is that the `.yaml` format is more human-readable.

An example command for running Peridigm:

````
Peridigm my_input.yaml
````

An example command for running Peridigm in parallel:

````
mpirun -np 4 Peridigm my_input.yaml
````

Peridigm is capable of operating on several types of discretizations, including Exodus/Genesis mesh files and user-generated text files. The most common use case is to generate a hexahedron and/or tetrahedron mesh using the [Cubit](https://cubit.sandia.gov/) mesh generator. In this case, hexahedron and/or tetrahedron elements are converted automatically by Peridigm into a meshless peridynamic discretization, i.e., a set of nodal volumes stored as (x, y, z, volume, block_id) data. Alternatively, users may supply a discretization via a text file in which (x, y, z, volume, block_id) data are stored directly. This option is provided for users who wish to construct discretizations, for example, using an in-house pre-processing code.

Parallel decomposition of Exodus/Genesis mesh files must be performed as a pre-processing step when running multi-processor Peridigm simulations. The decomp utility, provided by the SEACAS Trilinos package, provides a mechanism for decomposing an Exodus/Genesis mesh file.

An example command for partitioning an Exodus/Genesis mesh file with decomp:

````
decomp -p 4 my_mesh.g
````

Text file discretizations do not require this pre-processing step, they are partitioned automatically by Peridigm.

Peridigm generates output in the Exodus file format. The content of an Exodus output file is dictated by the Output section of a Peridigm input deck. Output may include primal quantities such a nodal displacements and velocities, as well as derived quantities such as stored elastic energy. The [ParaView](http://www.paraview.org/) visualization code is recommended for viewing Peridigm results. Additional options for parsing output data are available within the SEACAS Trilinos package.

The most effective way to learn how to use Peridigm is to run the example problems in the Peridigm/examples/ directory. These simulations were designed to highlight the most commonly-used features of Peridigm, including constitutive models, bond-failure rules, contact, explicit and implicit time integration, and I/O commands.

Questions regarding Peridigm should be sent to the [peridigm-users](https://software.sandia.gov/mailman/listinfo/peridigm-users) e-mail list.
