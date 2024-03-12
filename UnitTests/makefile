ARCH=${VLASIATOR_ARCH}
include ../../../MAKE/Makefile.${ARCH}

# These should be in Makefile.${ARCH}
# INC_GTEST=-I/projappl/project_2005018/lkotipal/libraries/gtest/include
# LIB_GTEST=-L/projappl/project_2005018/lkotipal/libraries/gtest/lib64 -lgtest -Wl,-rpath=/projappl/project_2005018/lkotipal/libraries/gtest/lib64

CXXFLAGS+=$(INC_MPI) $(INC_BOOST) $(INC_ZOLTAN) $(INC_GTEST) -DGTEST_HAS_MPI
LDFLAGS+=$(LIB_BOOST) $(LIB_ZOLTAN) $(LIB_GTEST)

all: simple.exe

simple.exe: simple.o
	$(CMP) $(CXXFLAGS) simple.o -o $@ $(LDFLAGS)

simple.o: simple.cpp common.hpp
	$(CMP) $(CXXFLAGS) -c simple.cpp

clean:
	rm *.exe
	rm *.o