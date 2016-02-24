# This Makefile exists just to cover basic usages and setup things we might want to do with this repository. For instance, installing required pip packages, running the tests, etc.

init:
	pip install -r requirements.txt
	pip install sphinx

# I'm using the syntax here for the pytest module, which integrates a bunch of useful testing ideas. 
# If we just use unittest or other variants, we can modify as needed.
test:
	py.test tests


docs:
	cd docs && make html

