

.PHONY: all
all: phantom

.PHONY: install
install:
	sudo python3 setup.py install

.PHONY: phantom
phantom: lib/pygfunc.pyf


.PHONY: dist
dist: phantom
	python3 setup.py sdist

.PHONY: upload
upload: dist
	twine upload --verbose dist/*

.PHONY: clean
clean:
	rm -f *.mod lib/*.o lib/pygfunc.pyf pygfunc.*.so
	rm -rf build dist
