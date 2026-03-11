Requirements
============

Binary Releases
---------------

Pre-compiled binaries of the latest release of Salmon for a number different
platforms are available available under the `Releases tab
<https://github.com/COMBINE-lab/salmon/releases>`_ of Salmon's `GitHub
repository <https://github.com/COMBINE-lab/salmon>`_.  You should be able to
get started quickly by finding a binary from the list that is compatible with
your platform.  Additionally, you can obtain a Docker image of the latest version
from DockerHub using:

::

    > docker pull combinelab/salmon
  

Requirements for Building from Source
-------------------------------------

* A C++17 conformant compiler
* CMake_ 3.24 or newer. Salmon uses the CMake build system to check, fetch and install
  dependencies, and to compile and install Salmon. CMake is available for all
  major platforms (though Salmon is currently unsupported on Windows.)
  
Installation
============

After downloading the Salmon source distribution and unpacking it, change into the top-level directory:

::

    > cd salmon

Salmon now provides ``CMakePresets.json`` presets for common workflows:

* ``dev`` for local development
* ``release`` for optimized builds
* ``asan`` for AddressSanitizer builds
* ``ci-system-deps`` for package-first CI
* ``ci-fetch-fallback`` for pinned fallback dependency builds

The preferred configure flow is:

::

    > cmake --preset dev


Salmon prefers system packages for dependencies and can fall back to pinned
source builds when ``SALMON_FETCH_MISSING_DEPS=ON``. Required backend
dependencies are:

* ``zlib-ng`` in compatibility mode for zlib functionality
* ``htslib`` for SAM/BAM/CRAM I/O
* ``mimalloc`` as the preferred allocator when ASan is disabled

Dependency resolution policy:

* ``SALMON_USE_SYSTEM_DEPS=ON`` (default): prefer ``find_package()`` and use
  installed packages when available.
* ``SALMON_FETCH_MISSING_DEPS=ON`` (default): fetch pinned fallbacks for
  missing required dependencies.
* ``SALMON_FETCH_MISSING_DEPS=OFF``: fail configure when required dependencies
  are not available from the system.

Supported CMake cache options are:

* ``SALMON_ENABLE_TESTS``
* ``SALMON_ENABLE_BENCHMARKS``
* ``SALMON_WARNINGS_AS_ERRORS``
* ``SALMON_ENABLE_LTO``
* ``SALMON_ENABLE_ASAN``
* ``SALMON_USE_SYSTEM_DEPS``
* ``SALMON_FETCH_MISSING_DEPS``
* ``SALMON_USE_ZLIB_NG``
* ``SALMON_USE_MIMALLOC``
* ``SALMON_USE_HTSLIB``

Legacy knobs like ``FETCH_BOOST`` and ``TBB_RECONFIGURE`` are not part of the
supported interface.

Examples:

development configure:

::

    > cmake --preset dev

release configure:

::

    > cmake --preset release

If you prefer to pass flags manually, you can still run the CMake configure
step as follows:

::
                                  
    > cmake [FLAGS] ..

The above command is the CMake configuration step. Next, run the build:

::

    > make

If the build is successful, the appropriate executables and libraries should be
created. Note that CMake has colored output by default, and link steps may
appear in red; this does not itself indicate a build failure.
                                  
Finally, after everything is built, the libraries and executable can be
installed with:

::
                                  
    > make install

You can test the installation by running

::

    > make test

This should run a simple test and tell you if it succeeded or not.

.. _CMake : http://www.cmake.org 
