# Use -Wno-missing-braces until GCC stops warning about valid C++ (http://gcc.gnu.org/bugzilla/show_bug.cgi?id=25137)
FLAGS = -O3 -W -Wall -Wextra -Wno-missing-braces -pedantic
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_mpi -lboost_program_options -lboost_serialization -lzoltan
CXX = mpic++ $(FLAGS)

all: simple

simple: simple.cpp ../../dccrg.hpp ../../dccrg_constant_geometry.hpp Makefile
	$(CXX) simple.cpp $(INCLUDES) -o simple

c: clean
clean:
	rm -f simple

