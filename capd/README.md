# CAPD

The implementation of the rigorous integration is found in this
directory. The implementation is given in
[`src/Q_zero.cpp`](/capd/src/Q_zero.cpp) and it compiles to an
executable `Q_zero` which is run by the Julia function `_Q_zero_capd`.

## Installation and building

To be able to run the computations the program first has to be
compiled. This requires installing the CAPD library and then compiling
the programs.

Note that the below instructions are for Linux. There is a lot of
variety between Linux systems, the description given below has worked
on several different systems but might need to be adapted.

### Installing CAPD

The library can be downloaded
[here](https://capd.sourceforge.net/capdRedHom/docs/html/capd_binary_install.html).
It is then installed using `configure` and `make`, as is common for C
and C++ libraries. The procedure that was used for the installation is
detailed below, which in this case installs it under `$HOME/capd`.

``` shell
# Create directory used for compilation
mkdir capd
cd capd

# Download and extract library
wget -O capd-5.3.0.tar.gz capd https://sourceforge.net/projects/capd/files/5.3.0/capd-5.3.0.tar.gz/download
tar xvf capd-5.3.0.tar.gz

# Configure and compile
mkdir build
cd build
../capd-5.3.0/configure --prefix=$HOME/capd --with-filib=internal --with-boost=internal --without-mpfr --without-gui
make -j
make install

# Output info about library and compiler
sha256sum ../capd-5.3.0.tar.gz
g++ --version
```

On the Agate HPC system (see [`HPC/README.d`](/HPC/README.md)) the
output from the last two lines were

``` shell
> sha256sum ../capd-5.3.0.tar.gz
e4100959a5409d330f8907d050f101a0485489075b4ce0d5eb2e349a2f8bf228  ../capd-5.3.0.tar.gz
> g++ --version
g++ (GCC) 8.2.0
Copyright (C) 2018 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
```

### Compiling programs

The programs `src/Q_zero.cpp` is compiled using the makefile found in
this directory. For the makefile to work it needs to have access to
the `capd-config` script that is installed with the CAPD library. The
script is found in the `bin/` folder of the CAPD installation
(`$HOME/capd/bin/` if the above installation procedure has been
followed). You can either make sure that `capd-config` is in the
`$PATH` or give the `CAPD_CONFIG` variable to make.

``` shell
# If capd-config is in the path you can just run make
make

# Otherwise give the path directly
CAPD_CONFIG=$HOME/capd/bin/capd-config make
```

You can see if the program is working by running

``` shell
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

It should give you the output

```
[-0.11224086758656471, -0.11224086758648889]
[-0.10823786762348438, -0.10823786762333835]
[-0.0075887592345378714, -0.0075887592339667275]
[0.018236639209484928, 0.01823663921004498]
```
