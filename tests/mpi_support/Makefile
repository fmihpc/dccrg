all:
	@$(MAKE) tests/mpi_support/executables -C ../..

t: test
test:
	@$(MAKE) tests/mpi_support/tests -C ../..

r: results
results:
	@$(MAKE) tests/mpi_support/results -C ../..

c: clean
clean:
	@$(MAKE) tests/mpi_support/clean -C ../..

