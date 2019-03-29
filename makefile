# If not running on gstar, then one can specify the path to the maps by calling: make MAP_PATH=/path/to/maps/without/quotes/and/ending/with/slash/
# To specify a different path you have to call "make clean" first.

# Prerequisites: CCfits, libpng, cuda, fftw
# The relevant lib and include directories must be in the corresponding environment path variables

CC    = g++
#CC_FLAGS = -std=c++11 -Wno-deprecated-gpu-targets # for a static library
CC_FLAGS   = -std=c++11 -fPIC
CC_LIBS    = -lpng -lCCfits -lcfitsio
INC   = -I include


CUDA_LOCAL = /usr/local/cuda-8.0/bin/nvcc
CUDA_LOCAL_PATH = /usr/local/cuda-8.0/lib64
CUDA_FLAGS = -std=c++11 --compiler-options '-fPIC' -Wno-deprecated-gpu-targets # for a dynamic library
CUDA_LIBS = -lcudart -lcufft

CUDA_OZSTAR = /apps/skylake/software/core/cuda/9.2.88/bin/nvcc
CUDA_OZSTAR_PATH = /apps/skylake/software/core/cuda/9.2.88/lib64


SOURCE_DIR = src
BUILD_DIR  = build
HEADER_DIR = include

SOURCES = image.cpp fixed_locs.cpp light_curve.cpp mpd.cpp profile.cpp velocities.cpp
MAG_MAP = magnification_map.cpp
GPU_SOURCES = gpu_functions.cu
CPU_SOURCES = cpu_functions.cpp


ifdef MAP_PATH
#QUOTED_MAP_PATH = $(addprefix "\",$(addsuffix \"",$(MAP_PATH)))
QUOTED_MAP_PATH = $(addprefix '",$(addsuffix "',$(MAP_PATH)))
MAP_PATH_FLAGS = -DMAP_PATH=$(QUOTED_MAP_PATH)
#$(info ${QUOTED_MAP_PATH})
endif
DOMAIN := $(shell dnsdomainname)





#SOURCES = $(shell find $(SOURCE_DIR) -type f -name '*.cpp')
#OBJECTS = $(patsubst $(SOURCE_DIR)/%,$(BUILD_DIR)/%,$(SOURCES:.cpp=.o))
#SOURCES += $(shell find $(SOURCE_DIR) -type f -name '*.cu')

PATH_SOURCES = $(patsubst %, $(SOURCE_DIR)/normal_source/%,$(SOURCES))
PATH_SOURCES += $(SOURCE_DIR)/magmap/magnification_map.cpp
TMP = $(patsubst $(SOURCE_DIR)/normal_source/%,$(BUILD_DIR)/%,$(PATH_SOURCES:.cpp=.o))
OBJECTS = $(patsubst $(SOURCE_DIR)/magmap/%,$(BUILD_DIR)/%,$(TMP:.cpp=.o))

PATH_GPU_SOURCES = $(patsubst %, $(SOURCE_DIR)/cpu_gpu/%,$(GPU_SOURCES))
GPU_OBJECTS = $(patsubst $(SOURCE_DIR)/cpu_gpu/%,$(BUILD_DIR)/%,$(PATH_GPU_SOURCES:.cu=.o))

PATH_CPU_SOURCES = $(patsubst %, $(SOURCE_DIR)/cpu_gpu/%,$(CPU_SOURCES))
CPU_OBJECTS = $(patsubst $(SOURCE_DIR)/cpu_gpu/%,$(BUILD_DIR)/%,$(PATH_CPU_SOURCES:.cpp=.o))

HEADERS = $(shell find $(HEADER_DIR) -type f -name '*.hpp')

#$(info $$PATH_SOURCES is [${PATH_SOURCES}])
#$(info $$PATH_GPU_SOURCES is [${PATH_GPU_SOURCES}])
#$(info $$PATH_CPU_SOURCES is [${PATH_CPU_SOURCES}])
#$(info $$OBJECTS is [${OBJECTS}])
#$(info $$GPU_OBJECTS is [${GPU_OBJECTS}])
#$(info $$CPU_OBJECTS is [${CPU_OBJECTS}])
#$(info $$HEADERS is [${HEADERS}])




# GPU PART
$(BUILD_DIR)/gpu_functions.o: $(SOURCE_DIR)/cpu_gpu/gpu_functions.cu $(HEADERS)
#	@echo $(DOMAIN
ifeq ($(DOMAIN),intra.astro.rug.nl)
	$(CUDA_LOCAL) $(CUDA_FLAGS) $(CUDA_LIBS) $(INC) -c -o $@ $< 
else
	$(CUDA_OZSTAR) $(CUDA_FLAGS) $(CUDA_LIBS) $(INC) -c -o $@ $< 
endif


# CPU PART
$(BUILD_DIR)/cpu_functions.o: $(SOURCE_DIR)/cpu_gpu/cpu_functions.cpp $(HEADERS)
	$(CC) $(CC_FLAGS) $(CC_LIBS) -lfftw3 $(INC) -c -o $@ $<

# MAP PATH
$(BUILD_DIR)/magnification_map.o: $(SOURCE_DIR)/magmap/magnification_map.cpp $(HEADERS)
	$(CC) $(CC_FLAGS) $(CC_LIBS) $(INC) $(MAP_PATH_FLAGS) -c -o $@ $<

# REST OF THE SOURCE CODE
$(BUILD_DIR)/%.o: $(SOURCE_DIR)/normal_source/%.cpp $(HEADERS)
	$(CC) $(CC_FLAGS) $(CC_LIBS) $(INC) -c -o $@ $< 



gpu: $(OBJECTS) $(GPU_OBJECTS) $(HEADERS)
ifeq ($(DOMAIN),intra.astro.rug.nl)
	g++ -shared -Wl,-soname,libgerlumph.so -L$(CUDA_LOCAL_PATH) -o lib/libgerlumph.so $(OBJECTS) $(GPU_OBJECTS) $(CC_LIBS) $(CUDA_LIBS)
else
	g++ -shared -Wl,-soname,libgerlumph.so -L$(CUDA_OZSTAR_PATH) -o lib/libgerlumph.so $(OBJECTS) $(GPU_OBJECTS) $(CC_LIBS) $(CUDA_LIBS) 
endif


cpu: $(OBJECTS) $(CPU_OBJECTS) $(HEADERS)
ifeq ($(DOMAIN),intra.astro.rug.nl)
	g++ -shared -Wl,-soname,libgerlumph.so -o lib/libgerlumph.so $(OBJECTS) $(CPU_OBJECTS) $(CC_LIBS) -lfftw3
else
	g++ -shared -Wl,-soname,libgerlumph.so -L/mnt/home/gvernard/myLibraries/fftw/lib -o lib/libgerlumph.so $(OBJECTS) $(CPU_OBJECTS) $(CC_LIBS) -lfftw3 
endif





test_gpu: gpu
	$(CC) -std=c++11 -I include -L lib -lgerlumph -o bin/readMap test/readMap.cpp


test_cpu: cpu
	$(CC) -std=c++11 -I include -L lib -lgerlumph -o bin/readMap test/readMap.cpp
#	$(CC) test/read_profile.cpp $(FLAGS) -I include -L lib $(LIBS) -o bin/read
#	$(CC) test/compare_profiles.cpp $(FLAGS) -I include -L lib $(LIBS) -o bin/compare
#	$(CC) test/main.cpp $(CC_FLAGS) -I include -L lib -lgerlumph -o bin/main
#	$(CC) test/create_profiles.cpp $(CC_FLAGS) -I include -L lib -lgerlumph -o bin/create


clean:	
	$(RM) -r $(BUILD_DIR)/* bin/* lib/*
