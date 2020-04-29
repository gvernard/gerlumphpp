#!/bin/bash

$PREFIX=/path/to/third-party-libraries/installation/directory

# Install fftw
wget http://www.fftw.org/fftw-3.3.8.tar.gz
tar -xvf fftw-3.3.8.tar.gz
mv fftw-3.3.8 fftw
cd fftw
./configure --prefix=$PREFIX/fftw --enable-shared
make CFLAGS=-fPIC
make install
cd ..

# Install cfitsio
wget http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-3.47.tar.gz
tar -xvf cfitsio-3.47.tar.gz
mv cfitsio-3.47 cfitsio
cd cfitsio
./configure --prefix=$PREFIX/cfitsio
make
make install
cd ..

# Install CCfits
wget https://heasarc.gsfc.nasa.gov/fitsio/CCfits/CCfits-2.5.tar.gz
tar -xvf CCfits-2.5.tar.gz
cd CCfits
./configure --prefix=$PREFIX/CCfits --with-cfitsio=$PREFIX/cfitsio
make
make install
cd ..

# Install libpng
wget https://sourceforge.net/projects/libpng/files/libpng16/1.6.37/libpng-1.6.37.tar.gz
tar -xvf libpng-1.6.37.tar.gz
cd libpng-1.6.37
./configure --prefix=$PREFIX/libpng
make check
sudo make install
cd ..



# OPTIONAL: Install jsoncpp
#########################################################################################################
git clone https://github.com/open-source-parsers/jsoncpp.git
cd jsoncpp
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=release -DBUILD_SHARED_LIBS=ON -DARCHIVE_INSTALL_DIR=. -DCMAKE_INSTALL_PREFIX=$PREFIX/jsoncpp -G "Unix Makefiles" .. 
make
sudo make install
cd ..



# ATTENTION!!! Add the following environment variables to your .bashrc file (or equivalent)
#########################################################################################################
CPATHT=$CPATH:$PREFIX/fftw/include:$PREFIX/cfitsio/include:$PREFIX/CCfits/include:$PREFIX/libpng/include:$PREFIX/jsoncpp/include
LIBRARY_PATH=$LIBRARY_PATH:$PREFIX/fftw/lib:$PREFIX/cfitsio/lib:$PREFIX/CCfits/lib:$PREFIX/libpng/lib:$PREFIX/jsoncpp/lib
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PREFIX/fftw/lib:$PREFIX/cfitsio/lib:$PREFIX/CCfits/lib:$PREFIX/libpng/lib:$PREFIX/jsoncpp/lib
