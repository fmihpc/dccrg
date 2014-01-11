# boost causes too many warnings at least on gcc-4.8.2 without -Wno-...
FLAGS = -O3 -W -Wall -Wextra -pedantic -std=c++0x -Wno-unused-local-typedefs
INCLUDES = -I$$HOME/include -L$$HOME/lib
CXX = mpic++ $(FLAGS)
RUN =

PROGRAMS = \
	enable_if \
	included \
	run_time

HEADERS = \
	../../dccrg.hpp \
	../../dccrg_get_cell_datatype.hpp

all: $(PROGRAMS)

enable_if: enable_if.cpp $(HEADERS) Makefile
	$(CXX) enable_if.cpp $(INCLUDES) -o enable_if

included: included.cpp $(HEADERS) Makefile
	$(CXX) included.cpp $(INCLUDES) -o included

run_time: run_time.cpp $(HEADERS) Makefile
	$(CXX) run_time.cpp $(INCLUDES) -o run_time

c: clean
clean:
	rm -f $(PROGRAMS)

t: test
test: enable_if run_time
	@for i in enable_if run_time; do echo -n "RUN $$i...  " && $(RUN) ./$$i && echo PASSED; done

