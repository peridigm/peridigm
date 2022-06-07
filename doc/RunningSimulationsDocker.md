# Running Simulations with Peridigm Docker Image

To run Peridigm from a [docker](https://docs.docker.com/) image, simply run the following command from the terminal command line of a computer with docker installed:

```bash
sudo docker run -v $PWD:/output peridigm/peridigm Peridigm [input_file.yaml]
```

where `[input_file.yaml]` should be replaced with the actual name of your
Peridigm input file.  This command should be run from the directory where `[input_file.yaml]` resides along with any geometry input files.

To launch a parallel Peridigm computation, we can use [docker-compose](https://docs.docker.com/compose/) (if installed) to first launch a network of containers running Peridigm with

```bash
sudo docker-compose up --scale peridigm=4 -d
```

from a directory that has a `docker-compose.yaml` file, e.g.
[examples/fragmenting_cylinder](examples/fragmenting_cylinder).  In this
example, we plan to use 4 processors, so we use the flag `--scale peridigm=4`.

Once the network is running, you can run Peridigm with


```bash
sudo docker-compose exec peridigm mpiexec -np 4 Peridigm [input_file.yaml]
```

where `[input_file.yaml]` should be replaced with the actual name of your
Peridigm input file.

This use of `sudo` in the above commands may or may not be needed depending on
you Docker installation, host operating system, and/or user configuration.
