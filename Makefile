
.PHONY: test

build:
	python setup.py build

develop:
	python setup.py develop

install:
	python setup.py install

test:
	coverage run -m py.test

clean:
	rm tests/data/*.fastq.?
