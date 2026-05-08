# CAPD

This directory contains the C++ implementation of the rigorous
integrator. The code is in [`src/Q_zero.cpp`](src/Q_zero.cpp), which
is compiled into the executable `Q_zero`. This executable is invoked
by the Julia function `_Q_zero_capd`.

## Installation and building

Before you can run the computations, you must compile the program.
This requires first installing the CAPD library and then building the
executable in this directory.

Note that the instructions below are for Linux. Since Linux systems
vary considerably, you may need to adapt them to your setup. However,
the procedure described below has worked on several different systems.

### Installing CAPD

You can download CAPD from the [CAPD GitHub
repository](https://github.com/CAPDGroup/CAPD), which also contains
installation instructions. The instructions below are for compiling
version 6.0.0. In this example, CAPD is installed under
`$HOME/CAPD`.

``` shell
# Download and extract library
wget -O CAPD-6.0.0.tar.gz https://github.com/CAPDGroup/CAPD/archive/refs/tags/v6.0.0.tar.gz
tar xvf CAPD-6.0.0.tar.gz
cd CAPD-6.0.0

# Configure and compile
mkdir build
cd build
# We disable multiprecision since it is not used.
cmake .. -DCAPD_ENABLE_MULTIPRECISION=false -DCMAKE_INSTALL_PREFIX=$HOME/CAPD
make -j
make install

# Output info about library and compiler
sha256sum ../../CAPD-6.0.0.tar.gz
g++ --version
```

On the Agate HPC system (see [`../HPC/README.md`](../HPC/README.md)),
the output from the last two lines was:

```shell
> sha256sum ../../CAPD-6.0.0.tar.gz
0ae254cb6477896c3c3f2b8cd1c4a6cb073fb1e1e08f71d0e0bd2587d60c9620  ../../CAPD-6.0.0.tar.gz
> g++ --version
g++ (GCC) 8.5.0 20210514 (Red Hat 8.5.0-28)
Copyright (C) 2018 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
```

### Building the executable

The file `src/Q_zero.cpp` is compiled using the `Makefile` in this
directory. For the build to work, `make` must be able to find the
`capd-config` script installed with the CAPD library. This script is
located in the `bin/` directory of the CAPD installation
(`$HOME/CAPD/bin/` if you followed the installation procedure above).
You can either make sure that `capd-config` is in your `$PATH` or pass
the `CAPD_CONFIG` variable to `make`.

```shell
# If capd-config is in your PATH, you can just run make
make

# Otherwise, give the path directly
CAPD_CONFIG=$HOME/CAPD/bin/capd-config make
```

You can verify that the program is working by running:

```shell
echo '[1.232037525342875, 1.232037525342875]
[-0.0, 0.0]
[-0.0, 0.0]
[-0.0, 0.0]
1
[0.8531088349377206, 0.8531088349377206]
[-0.0, 0.0]
[1.0, 1.0]
[2.3, 2.3000000000000003]
[-0.0, 0.0]
[-0.0, 0.0]
[10.0, 10.0]
0
0
0
1.0e-11
' | ./build/Q_zero
```

It should produce the following output:

```
[-0.11224086758656471, -0.11224086758648889]
[-0.10823786762348438, -0.10823786762333835]
[-0.0075887592345378714, -0.0075887592339667275]
[0.018236639209484928, 0.01823663921004498]
```
