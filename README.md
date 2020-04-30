# gerlumphpp

This is a cpp library of functions to perform various microlensing tasks, e.g. read maps in a given format, doing convolutions, extracting light curves, etc.
Optionally, the convolutions can be performed on the GPU using cuda.



## Getting Started
### Prerequisites

There is a number of third-party libraries required to install gerlumphpp, namely: **fftw3**, **cfitsio**, **CCfits**, and **libpng** (optionally **jsoncpp**).
These are required for the CPU-based version.
The GPU requires additionally **cuda** (nvcc compiler) and **cufft**.
For the latter, a root -based system-wide installation is preferred, while the rest can be installed locally.
The **prerequisites.sh** script takes care of a local installation by replacing the variable **$PREFIX** with your own local path and typing:

```
bash prerequisites.sh
```

This script will download and install the CPU-based third-party software required.
It only requires **cmake** to be already installed -  apretty standard package on Unix-like systems.
Only two *make install* commands in this script require root permission (sudo), which, however, can be avoided by changing the write permissions on the **$PREFIX** directory.

As a final step, the environment variables **CPATH**, **LIBRARY_PATH**, and **LD_LIBRARY_PATH** at the end of the script need to be included somewhere in the user's path, e.g. in a .bashrc file, or equivalent.



### Installing

To install gerlumphpp one **must** specify an absolute path to magnification maps by passing the **MAP_PATH** argument to make.
Then, specify which version of the library to build by typing **cpu** or **gpu**.
An example make is (at the root directory of gerlumphpp):

```
make MAP_PATH=/path/to/maps/without/quotes/and/ending/with/slash/ cpu
```

To specify a different **MAP_PATH** one has to call "make clean" first.

Finally, to be able to compile applications using gerlumphpp without specifying the path to *include* and *lib*, add the following in your .bashrc (similarly for other shells):

```
CPATH=$CPATH:/path/to/gerlumphpp/include
LIBRARY_PATH=$LIBRARY_PATH:/path/to/gerlumphpp/lib
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/gerlumphpp/lib
```





## Running the tests

In tests/example.cpp there is a complete example of reading a map, convolving with extended source light profiles, and extracting light curves or magnification probability distributions.
To compile the tests type:

```
make tests_cpu
```
or
```
make tests_gpu
```
according to which version of the library you have compiled.

To run the example, [download]() a map from GERLUMPH and place it at your specified **MAP_PATH**.
You will need a directory named "12345" (for this specific example, but it can be whatever) that contains the *map.bin* and *mapmeta.dat* files.
Then type:

```
cd tests
./example
```

and the example code will produce a list of output images and data.





## Authors

**Georgios Vernardos** ( [github](https://github.com/gvernard)  - [homepage](http://astronomy.swin.edu.au/~gvernard/) )



