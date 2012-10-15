========
PyLSMLIB
========

PyLSMLIB_ is a Python package that wraps the LSMLIB_ level set library
using Cython_. It arose out of the need to replace the rudimentary
level set capability in FiPy_ with a faster C-based level set library.
The interface is designed to be pythonic and match the interface to
Scikit-fmm_.

Hosted
======

<https://github.com/ktchu/LSMLIB.git/pylsmlib>

PyLSMLIB_ is hosted at Github_ in a subdirectory of LSMLIB_.

Requirements
============

The main requirements are working versions of LSMLIB_ and
Numpy_. Cython_ should not be a requirement as the generated ``*.c``
files are included with the distribution. See `requirements.txt`_ for
specific versions of dependencies used to run this package on the
maintainer's system. The `requirements.txt`_ file is auto-generated so
most of the packages listed are not necessarily required. However, if
you are having issues with installation and missing packages this
would be a good place to start looking.

Installation
============

The following should get you most of the way to an install.

::

    $ pip install numpy

Clone the git repository.

::

    $ git clone git://github.com/ktchu/LSMLIB.git LSMLIB

See the Github_ project page for further details. After cloning,
install LSMLIB_ (consult `LSMLIB's install`_ if you have issues).

::
  
    $ cd .../LSMLIB
    $ mkdir build
    $ cd build
    $ ../configure
    $ make
    $ make install

To install PyLSMLIB_.

::

    $ cd .../LSMLIB/pylsmlib
    $ python setup.py install

Testing
=======

To run the tests

::

    >>> import pylsmlib
    >>> pylsmlib.test()

Documentation
=============

To generate the PyLSMLIB_ documentation.

::

    $ cd .../LSMLIB/pylsmlib/doc
    $ make html

.. _LSMLIB: http://ktchu.serendipityresearch.org/software/lsmlib/index.html
.. _PyLSMLIB: https://github.com/ktchu/LSMLIB/tree/master/pylsmlib
.. _Github: https://github.com/ktchu/LSMLIB
.. _requirements.txt: https://github.com/ktchu/LSMLIB/blob/master/pylsmlib/requirements.txt
.. _Cython: http://cython.org/
.. _FiPy: http://www.ctcms.nist.gov/fipy/
.. _Scikit-fmm: http://packages.python.org/scikit-fmm/
.. _Numpy: http://numpy.scipy.org/
.. _LSMLIB's install: https://github.com/ktchu/LSMLIB/blob/master/INSTALL
