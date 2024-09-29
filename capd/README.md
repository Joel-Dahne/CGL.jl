# Installing CAPD
All below instructions are for Linux.

The library can be downloaded
[here](https://capd.sourceforge.net/capdRedHom/docs/html/capd_binary_install.html).
It is then installed using `configure` and `make` in the standard way.
The below code has worked for me on several Linux systems.

``` shell
# Download and extract library
wget https://sourceforge.net/projects/capd/files/5.3.0/capd-5.3.0.tar.gz/download
tar xvf capd-5.3.0.tar.gz

# Configure and compile
cd capd-5.3.0
./configure --prefix=$HOME/capd --with-filib=no --with-mpfr=no --with-gui=no
make -j
make install
```

For running the code make sure that the CAPD library is in your
`LD_LIBRARY_PATH`. In my case that is `$HOME/Programs/capd/lib/`.

For Agate the compilation was done using GCC 11.3.0.
