CUDA   = nvcc
FLAGS = -std=c++11

ODIR = obj
LIBS = -lcufft -lpng -lCCfits
DEPS       = magnification_map.hpp mpd.hpp profile.hpp kernel.hpp light_curve.hpp fixed_locs.hpp image.hpp

OBJ_MAIN   = magnification_map.o   mpd.o   profile.o   kernel.o   light_curve.o   fixed_locs.o   image.o   main.o
OBJ_CREATE = magnification_map.o   mpd.o   profile.o   kernel.o   light_curve.o   fixed_locs.o   image.o   create_profile.o
OBJ_READ   = magnification_map.o   mpd.o   profile.o   kernel.o   light_curve.o   fixed_locs.o   image.o   read_profile.o

#Pad object names with object dir
_OBJ_MAIN = $(patsubst %,$(ODIR)/%,$(OBJ_MAIN))
_OBJ_CREATE = $(patsubst %,$(ODIR)/%,$(OBJ_CREATE))
_OBJ_READ = $(patsubst %,$(ODIR)/%,$(OBJ_READ))


$(ODIR)/%.o: %.cpp $(DEPS)
	$(CUDA) -c -o $@ $< $(FLAGS)

$(ODIR)/%.o: %.cu $(DEPS)
	$(CUDA) -c -o $@ $< $(FLAGS)


all: main create read

main: $(_OBJ_MAIN)
	$(CUDA) -o $@ $^ $(FLAGS) $(LIBS)

create: $(_OBJ_CREATE)
	$(CUDA) -o $@ $^ $(FLAGS) $(LIBS)

read: $(_OBJ_READ)
	$(CUDA) -o $@ $^ $(FLAGS) $(LIBS)


