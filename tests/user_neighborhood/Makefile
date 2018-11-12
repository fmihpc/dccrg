all:
	@$(MAKE) tests/user_neighborhood/executables -C ../..

t: test
test:
	@$(MAKE) tests/user_neighborhood/tests -C ../..

r: results
results:
	@$(MAKE) tests/user_neighborhood/results -C ../..

c: clean
clean: results
	@$(MAKE) tests/user_neighborhood/clean -C ../..

