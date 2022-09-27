# gerlumphpp

This is a C++ library of functions to perform various microlensing tasks, e.g. read maps, perform convolutions, extract light curves, etc.
Optionally, the convolutions can be performed on the GPU using cuda, which leads to a approx. 5 times speedup.



## Getting Started
### Prerequisites

There is a number of third-party libraries required to install gerlumphpp, namely: **fftw3**, **cfitsio**, **CCfits**, and **libpng**.
These are required for the CPU-based version.
To enable GPU support **cuda** (the nvcc compiler, tested with 11.7, does not work with 11.5 due to a bug) and **cufft** are required as well.
For the latter two, a root-based system-wide installation is preferred, while the rest can be installed locally.
Finally, **autotools** are required to perform the installation (the autoreconf command). 


### Installing

To install gerlumphpp one needs to call `autoreconf -i' and then the usual './configure, make, make install'.
However, there are some mandatory options that need to be passed to the configure script.

- --with-map-path: the user must specify the absolute path to a directory containing GERLUMPH maps
- --enable-gpu: either 'yes' or 'no', that will compile the library with or without GPU support

In addtion, if the required third-party libraries listed above are not installed in a standard system location, they will need to be explicitly provided.
This is possible via '--with-<library_name>=/path/to/library' options passed to the configure script.
An example call to the configure script could look like the following:

```
./configure --prefix=/path/to/installation/of/gerlumphpp --with-fftw3=/path/to/libraries/fftw --with-cfitsio=/path/to/libraries/cfitsio --with-CCfits=/path/to/libraries/CCfits --with-png=/path/to/libraries/libpng --enable-gpu=yes --with-map-path=/path/to/gerlumph/maps/
```


### Enabling/disabling GPU support

One may need to re-compile gerlumphpp at some point because the location of the third-party libraries or of the GERLUMPH maps has changed.
Most commonly though, one may want to enable/disable GPU support depending on the system configuration.
In this case, the same steps as above need to be followed with the appropriate options set.
In fact, one can have two versions of the library, one for the CPU and one for the GPU, simultaneously present in the system installed at different locations. 


### Finalizing

For convenience, in order to be able to compile programs that use gerlumphpp without the need to explicitly specify the path to the headers and the libarry (e.g. the -I, -L, and -Wl,-rpath options to the g++ compiler and linker), one can define the following environment variables:

```
CPATH=$CPATH:/path/to/installation/of/gerlumphpp/include
LIBRARY_PATH=$LIBRARY_PATH:/path/to/installation/of/gerlumphpp/lib
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/installation/of/gerlumphpp/lib
```

For the bash shell, one can export these variables from within the .bashrc file.






## Running the tests

In tests/example.cpp there is a complete example of reading a map, convolving with extended source light profiles, and extracting light curves or magnification probability distributions.
To run the example, [download]() a map from GERLUMPH and place it at your specified **MAP_PATH**.
You will need a directory named "12345" (for this specific example, but it can be whatever) that contains the *map.bin* and *mapmeta.dat* files.
Then compile the example linking it properly with the gerlumphpp library and run it.
The example program will produce a list of output images and additional data like light curves.



## Authors

**Georgios Vernardos** ( [github](https://github.com/gvernard)  - [homepage](http://astronomy.swin.edu.au/~gvernard/) )



