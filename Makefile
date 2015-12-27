# to compile and test pick one makefile from makefiles directory,
# or write a new one for your environment, and put it here:
ENVIRONMENT_MAKEFILE = makefiles/homedir

#
# The lines below are not intended to be modified by users
#
CXXFLAGS = -std=c++11 -W -Wall -Wextra -pedantic -O3
include $(ENVIRONMENT_MAKEFILE)

# filled by project makefiles
EXECUTABLES =
TESTS =
RESULTS =
CLEAN =

default: all

DCCRG_HEADERS = \
  dccrg_cartesian_geometry.hpp \
  dccrg_get_cell_datatype.hpp \
  dccrg.hpp \
  dccrg_iterator_support.hpp \
  dccrg_length.hpp \
  dccrg_mapping.hpp \
  dccrg_mpi_support.hpp \
  dccrg_no_geometry.hpp \
  dccrg_stretched_cartesian_geometry.hpp \
  dccrg_topology.hpp \
  dccrg_types.hpp

include \
  examples/project_makefile \
  tests/init/project_makefile \
  tests/get_cell_datatype/project_makefile


all: $(EXECUTABLES)

t: test
test: all $(TESTS)

# removes simulation results
r: results
results: $(RESULTS)

c: clean
clean: results $(CLEAN)


# Rules to run tests common to all projects
%.tst: %.exe
	@printf RUN\ $<...\ \  && $(RUN) ./$< && printf "PASS\n" && touch $@

%.mtst: %.exe
	@printf MPIRUN\ $<...\ \  && $(MPIRUN) ./$< && printf "PASS\n" && touch $@

