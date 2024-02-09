#!/bin/bash
LIB_PATH=/home/giorgos/myLibraries/libraries

g++ -o example example.cpp \
    -I${LIB_PATH}/gerlumphpp/include/ \
    -I${LIB_PATH}/CCfits/include/ \
    -I${LIB_PATH}/cfitsio/include/ \
    -L${LIB_PATH}/gerlumphpp/lib -Wl,-rpath,${LIB_PATH}/gerlumphpp/lib/ \
    -L${LIB_PATH}/CCfits/lib -Wl,-rpath,${LIB_PATH}/CCfits/lib/ \
    -L${LIB_PATH}/cfitsio/lib -Wl,-rpath,${LIB_PATH}/cfitsio/lib/ \
    -lgerlumph
