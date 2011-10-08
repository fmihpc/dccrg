# Use -Wno-missing-braces until GCC stops warning about valid C++ (http://gcc.gnu.org/bugzilla/show_bug.cgi?id=25137)
FLAGS = -O3 -W -Wall -Wextra -pedantic -Wno-missing-braces -std=c++0x
#FLAGS = -O0 -g -W -Wall -Wextra -pedantic -std=c++0x
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_mpi -lboost_serialization -lzoltan
CXX = mpic++ $(FLAGS)

all: load_balancing_test hierarchical

load_balancing_test: load_balancing_test.cpp ../../dccrg.hpp ../../dccrg_constant_geometry.hpp Makefile
	$(CXX) load_balancing_test.cpp $(INCLUDES) -o load_balancing_test

hierarchical: hierarchical.cpp ../../dccrg.hpp ../../dccrg_constant_geometry.hpp Makefile
	$(CXX) hierarchical.cpp $(INCLUDES) -o hierarchical

c: clean
clean:
	rm -f load_balancing_test hierarchical *vtk *visit
