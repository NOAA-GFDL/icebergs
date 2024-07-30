Building the documentation
==========================

Using just doxygen
```
doxygen
```
will create HTML in the directory `html/`. Error messages will appear in `doxygen.log`.

The readthedocs site, https://kid.readthedocs.io/, is generated using sphinx (see https://sphinx-doc.org/, install with `conda install sphinx`) with
```
make html
```
