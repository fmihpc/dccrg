ARCH=${VLASIATOR_ARCH}
include ../../../MAKE/Makefile.${ARCH}

# These should be in Makefile.${ARCH}
# INC_GTEST=-I/projappl/project_2005018/lkotipal/libraries/gtest/include
# LIB_GTEST=-L/projappl/project_2005018/lkotipal/libraries/gtest/lib64 -lgtest -Wl,-rpath=/projappl/project_2005018/lkotipal/libraries/gtest/lib64
# LIB_NOPROFILE = -L$(LIBRARY_PREFIX)/phiprof/lib -lnophiprof -Wl,-rpath=$(LIBRARY_PREFIX)/phiprof/lib

CXXFLAGS+=$(INC_MPI) $(INC_BOOST) $(INC_ZOLTAN) $(INC_PROFILE) $(INC_GTEST) -DGTEST_HAS_MPI
LDFLAGS+=$(LIB_BOOST) $(LIB_ZOLTAN) $(LIB_GTEST) $(LIB_NOPROFILE)

TESTS=validation.exe grid.exe

all: $(TESTS)

%.exe: %.o
	$(CMP) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CMP) -c $(CXXFLAGS) $< -o $@

clean:
	rm $(TESTS)
