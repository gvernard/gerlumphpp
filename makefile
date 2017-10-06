CUDA  = nvcc
FLAGS = -std=c++11 -Wno-deprecated-gpu-targets
LIBS  = -lcufft -lpng -lCCfits -lgerlumph
INC   = -I include

SOURCE_DIR = src
BUILD_DIR  = build
HEADER_DIR = include


SOURCES = $(shell find $(SOURCE_DIR) -type f -name *.cpp)
SOURCES += $(shell find $(SOURCE_DIR) -type f -name *.cu)
TMP = $(patsubst $(SOURCE_DIR)/%,$(BUILD_DIR)/%,$(SOURCES:.cpp=.o))
OBJECTS = $(TMP:.cu=.o)
HEADERS = $(shell find $(HEADER_DIR) -type f -name *.hpp)


#$(info $$SOURCES is [${SOURCES}])
#$(info $$OBJECTS is [${HEADERS}])




libgerlumph.a: $(OBJECTS) $(HEADERS)
	ar rcs lib/$@ $(OBJECTS)


build/magnification_map.o: src/magnification_map.cu $(HEADERS)
	$(CUDA) $(FLAGS) $(INC) -c -o $@ $< 

build/%.o: src/%.cpp $(HEADERS)
	$(CUDA) $(FLAGS) $(INC) -c -o $@ $< 


test: libgerlumph.a
	$(CUDA) test/create_profiles.cpp $(FLAGS) -I include -L lib $(LIBS) -o bin/create
#	$(CUDA) test/read_profile.cpp $(FLAGS) -I include -L lib $(LIBS) -o bin/read
#	$(CUDA) test/compare_profiles.cpp $(FLAGS) -I include -L lib $(LIBS) -o bin/compare
	$(CUDA) test/main.cpp $(FLAGS) -I include -L lib $(LIBS) -o bin/main


clean:	
	$(RM) -r $(BUILD_DIR)/* bin/* lib/*
