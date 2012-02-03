# Use -Wno-missing-braces until GCC stops warning about valid C++ (http://gcc.gnu.org/bugzilla/show_bug.cgi?id=25137)
FLAGS = -O3 -W -Wall -Wextra -pedantic -Wno-missing-braces
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_mpi -lboost_program_options -lboost_serialization -lzoltan
CXX = mpic++ $(FLAGS)

PROGRAMS = 2d 2d_debug dc2vtk

HEADERS = \
	../../dccrg.hpp \
	../../dccrg_index.hpp \
	../../dccrg_constant_geometry.hpp

all: $(PROGRAMS)

2d: 2d.cpp $(HEADERS) Makefile
	$(CXX) -DNDEBUG 2d.cpp $(INCLUDES) -o 2d

2d_debug: 2d.cpp $(HEADERS) Makefile
	$(CXX) -DDEBUG 2d.cpp $(INCLUDES) -o 2d_debug

dc2vtk: dc2vtk.cpp ../../dccrg_constant_geometry.hpp Makefile
	$(CXX) dc2vtk.cpp $(INCLUDES) -o dc2vtk

d: data
data:
	rm -f *.dc *.vtk

c: clean
clean: data
	rm -f $(PROGRAMS)

