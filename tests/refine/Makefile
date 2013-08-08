FLAGS = -O3 -W -Wall -Wextra -pedantic -std=c++0x
INCLUDES = -I../.. -I$$HOME/include -L$$HOME/lib -lboost_mpi -lboost_program_options -lboost_serialization -lzoltan
CXX = mpic++ $(FLAGS)

all: refine_simple unrefine_simple scalability

refine_simple: refine_simple.cpp ../../dccrg.hpp Makefile
	$(CXX) -DDCCRG_TRANSFER_USING_BOOST_MPI refine_simple.cpp $(INCLUDES) -o refine_simple

unrefine_simple: unrefine_simple.cpp ../../dccrg.hpp Makefile
	$(CXX) -DDCCRG_TRANSFER_USING_BOOST_MPI unrefine_simple.cpp $(INCLUDES) -o unrefine_simple

scalability: scalability.cpp ../../dccrg.hpp Makefile
	$(CXX) -DDCCRG_TRANSFER_USING_BOOST_MPI scalability.cpp $(INCLUDES) -o scalability

c: clean
clean:
	rm -f refine_simple unrefine_simple scalability *.visit *.vtk
