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

* A C++11 conformant compiler (currently tested with GCC>=4.7 and Clang>=3.4)
* CMake_. Salmon uses the CMake build system to check, fetch and install
  dependencies, and to compile and install Salmon. CMake is available for all
  major platforms (though Salmon is currently unsupported on Windows.)
  
Installation
============

After downloading the Salmon source distribution and unpacking it, change into the top-level directory:

::

    > cd salmon

Then, create and out-of-source build directory and change into it:

::

    > mkdir build
    > cd build


Salmon makes extensive use of Boost_.  We recommend installing the most
recent version (1.55) systemwide if possible. If Boost is not installed on your
system, the build process will fetch, compile and install it locally.  However,
if you already have a recent version of Boost available on your system, it make
sense to tell the build system to use that.

If you have Boost installed you can tell CMake where to look for it. Likewise,
if you already have `Intel's Threading Building Blocks
<http://threadingbuildingblocks.org/>`_ library installed, you can tell CMake
where it is as well. The flags for CMake are as follows:

* -DFETCH_BOOST=TRUE --  If you don't have Boost installed (or have an older
  version of it), you can provide the FETCH_BOOST flag instead of the
  BOOST_ROOT variable, which will cause CMake to fetch and build Boost locally.

* -DBOOST_ROOT=<boostdir> -- Tells CMake where an existing installtion of Boost
  resides, and looks for the appropritate version in <boostdir>.  This is the
  top-level directory where Boost is installed (e.g. /opt/local).

* -DTBB_INSTALL_DIR=<tbbroot> -- Tells CMake where an existing installation of
  Intel's TBB is installed (<tbbroot>), and looks for the apropriate headers
  and libraries there. This is the top-level directory where TBB is installed
  (e.g. /opt/local).

* -DCMAKE_INSTALL_PREFIX=<install_dir> -- <install_dir> is the directory to
  which you wish Salmon to be installed.  If you don't specify this option,
  it will be installed locally in the top-level directory (i.e. the directory
  directly above "build").

There are a number of other libraries upon which Salmon depends, but CMake 
should fetch these for you automatically.

Setting the appropriate flags, you can then run the CMake configure step as
follows:

::
                                  
    > cmake [FLAGS] ..

The above command is the cmake configuration step, which *should* complain if
anything goes wrong.  Next, you have to run the build step. Depending on what
libraries need to be fetched and installed, this could take a while
(specifically if the installation needs to install Boost).  To start the build,
just run make.

::

    > make

If the build is successful, the appropriate executables and libraries should be
created. There are two points to note about the build process.  First, if the
build system is downloading and compiling boost, you may see a large number of
warnings during compilation; these are normal.  Second, note that CMake has
colored output by default, and the steps which create or link libraries are
printed in red.  This is the color chosen by CMake for linking messages, and
does not denote an error in the build process. 
                                  
Finally, after everything is built, the libraries and executable can be
installed with:

::
                                  
    > make install

You can test the installation by running

::

    > make test

This should run a simple test and tell you if it succeeded or not.

.. _CMake : http://www.cmake.org 
.. _Boost: http://www.boost.org
