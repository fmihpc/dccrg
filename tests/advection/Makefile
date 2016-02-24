all:
	@$(MAKE) tests/advection/executables -C ../..

t: test
test:
	@$(MAKE) tests/advection/tests -C ../..

r: results
results:
	@$(MAKE) tests/advection/results -C ../..

c: clean
clean: results
	@$(MAKE) tests/advection/clean -C ../..

