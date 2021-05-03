# Directory overview

| Directory     | Description                                         |
| --------------| --------------------------------------------------- |
| test_2d       | Directory for building 2D EB executable             |
| test_3d       | Directory for building 3D EB executable             |
| test_no_eb_2d | Directory for building 2D non-EB executable         |
| test_no_eb    | Directory for building 3D non-EB executable         |
| src           | Source files                                        |


# Using incflo

To build and run incflo, you will need to fist clone the AMReX and AMReX-Hydro repositories.

We suggest using the _main_ branch of both repositories.

## Check out AMReX 

Clone AMReX from the official Git repository
```shell
> git clone git@github.com:AMReX-Codes/amrex.git
```

Clone AMReX-Hydro from the official Git repository
```shell
> git clone git@github.com:AMReX-Codes/AMReX-Hydro.git
```

## Build and run an example incflo problem
Clone and build incflo
```shell
> git clone git@github.com:AMReX-Codes/incflo.git
> cd incflo/test_3d
> make -j4
> mpirun -np 4 incflo3d.gnu.MPI.ex benchmark.channel_cylinder-x
```

# Contributing

We welcome contributions in the form of pull-requests from anyone.
