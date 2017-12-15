# If not running on gstar, then one can specify the path to the maps by calling: make MAP_PATH=/path/to/maps/without/quotes/and/ending/with/slash/
# To specify a different path you have to call "make clean" first.

CC    = g++
CUDA  = nvcc
#CC_FLAGS = -std=c++11 -Wno-deprecated-gpu-targets # for a static library
CC_FLAGS   = -std=c++11 -fPIC
CUDA_FLAGS = -std=c++11 --compiler-options '-fPIC' -Wno-deprecated-gpu-targets # for a dynamic library

ifdef MAP_PATH
QUOTED_MAP_PATH = $(addprefix "\",$(addsuffix \"",$(MAP_PATH)))
CC_FLAGS += -DMAP_PATH=$(QUOTED_MAP_PATH)
CUDA_FLAGS += -DMAP_PATH=$(QUOTED_MAP_PATH)
#$(info ${QUOTED_MAP_PATH})
endif


CC_LIBS    = -lpng -lCCfits
CUDA_LIBS  = -lcuda -lcudart -lcufft
INC   = -I include

SOURCE_DIR = src
BUILD_DIR  = build
HEADER_DIR = include


SOURCES = $(shell find $(SOURCE_DIR) -type f -name *.cpp)
#OBJECTS = $(patsubst $(SOURCE_DIR)/%,$(BUILD_DIR)/%,$(SOURCES:.cpp=.o))
SOURCES += $(shell find $(SOURCE_DIR) -type f -name *.cu)
TMP = $(patsubst $(SOURCE_DIR)/%,$(BUILD_DIR)/%,$(SOURCES:.cpp=.o))
OBJECTS = $(TMP:.cu=.o)
HEADERS = $(shell find $(HEADER_DIR) -type f -name *.hpp)


#$(info $$SOURCES is [${SOURCES}])
#$(info $$OBJECTS is [${HEADERS}])



DOMAIN := $(shell dnsdomainname)
libgerlumph.so: $(OBJECTS) $(HEADERS)
ifeq ($(DOMAIN),hpc.swin.edu.au)
	g++ -shared -Wl,-soname,libgerlumph.so -L/usr/local/cuda-7.5/lib64 $(CC_LIBS) $(CUDA_LIBS) -o lib/libgerlumph.so $(OBJECTS)
else
	g++ -shared -Wl,-soname,libgerlumph.so -L/usr/local/cuda/lib64 $(CC_LIBS) $(CUDA_LIBS) -o lib/libgerlumph.so $(OBJECTS)
endif


build/magnification_map.o: src/magnification_map.cu $(HEADERS)
	$(CUDA) $(CUDA_FLAGS) $(CUDA_LIBS) $(INC) -c -o $@ $< 

build/%.o: src/%.cpp $(HEADERS)
	$(CC) $(CC_FLAGS) $(CC_LIBS) $(INC) -c -o $@ $< 


test: libgerlumph.so
#	$(CC) test/read_profile.cpp $(FLAGS) -I include -L lib $(LIBS) -o bin/read
#	$(CC) test/compare_profiles.cpp $(FLAGS) -I include -L lib $(LIBS) -o bin/compare
	$(CC) test/main.cpp $(CC_FLAGS) -I include -L lib -lgerlumph -o bin/main
	$(CC) test/create_profiles.cpp $(CC_FLAGS) -I include -L lib -lgerlumph -o bin/create


clean:	
	$(RM) -r $(BUILD_DIR)/* bin/* lib/*
