

.PHONY: all
all: phantom

.PHONY: install
install:
	sudo python2 setup.py install

.PHONY: phantom
phantom: lib/pygfunc.c
	python setup.py build_ext --inplace

lib/pygfunc.c:
	cython -a lib/pygfunc.pyx

.PHONY: dist
dist:
	python2 setup.py sdist

.PHONY: upload
upload: dist
	twine upload --verbose dist/*

.PHONY: clean
clean:
	rm -f *.mod lib/*.o lib/pygfunc.c lib/pygfunc.html pygfunc.so
	rm -rf build dist
