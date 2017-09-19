CUDA   = nvcc
FLAGS = -std=c++11

ODIR = obj
LIBS = -lcufft -lpng
DEPS = magnification_map.hpp mpd.hpp profile.hpp kernel.hpp light_curve.hpp
OBJ  = magnification_map.o   mpd.o   profile.o   kernel.o   light_curve.o   main.o

#Pad object names with object dir
_OBJ = $(patsubst %,$(ODIR)/%,$(OBJ))


$(ODIR)/%.o: %.cpp $(DEPS)
	$(CUDA) -c -o $@ $< $(FLAGS)

$(ODIR)/%.o: %.cu $(DEPS)
	$(CUDA) -c -o $@ $< $(FLAGS)

main: $(_OBJ)
	$(CUDA) -o $@ $^ $(FLAGS) $(LIBS)


