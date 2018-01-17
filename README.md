# Open Porous Media Core Library

THIS MODULE IS DEPRECATED.
The code is now integrated in the other modules.

These are release notes for opm-core.


CONTENT
-------

opm-core is the core library within OPM and contains the following 

* Fluid properties (basic PVT models and rock properties)
* Grid handling (cornerpoint grids, unstructured grid interface)
* Linear Algebra (interface to different linear solvers)
* Pressure solvers (various discretization schemes, flow models)
* Simulators (some basic examples of simulators based on sequential splitting schemes)
* Transport solvers (various discretization schemes, flow models)
* Flow diagnostics (time-of-flight and tracer solvers, diagnostic functions)
* Utilities (input and output processing, unit conversion)
* Wells (basic well handling)


LICENSE
-------

The library is distributed under the GNU General Public License,
version 3 or later (GPLv3+).


PLATFORMS
---------

The opm-core module is designed to run on Linux platforms. It is also
regularly run on Mac OS X. No efforts have been made to ensure that
the code will compile and run on windows platforms.


DEPENDENCIES FOR DEBIAN BASED DISTRIBUTIONS (Debian Squeeze/Ubuntu Precise)
---------------------------------------------------------------------------

# packages necessary for building
sudo apt-get install -y build-essential gfortran cmake cmake-data util-linux

# packages necessary for documentation
sudo apt-get install -y doxygen ghostscript texlive-latex-recommended pgf

# packages necessary for version control
sudo apt-get install -y git-core

# basic libraries necessary for both DUNE and OPM
sudo apt-get install -y libboost-all-dev libsuperlu3-dev libsuitesparse-dev

# for server edition of Ubuntu add-apt-repository depends on
sudo apt-get install python-software-properties

# add this repository for necessary backports (required for Ubuntu Precise)
sudo add-apt-repository -y ppa:opm/ppa
sudo apt-get update

# parts of DUNE needed
sudo apt-get install libdune-common-dev libdune-istl-dev libdune-grid-dev

# libraries necessary for OPM
sudo apt-get install -y libtinyxml-dev

# Ensemble based Reservoir Tool Eclipse utilities module
# IMPORTANT: if you install this (binary) version of ERT,
#            you will get the 2015.04 release version. That
#            is only compatible with the 2015.04 release version
#            of OPM! If you are building OPM from source you should
#            use the latest master branches of both ERT and OPM.
sudo apt-get install ert.ecl

Note: You should compile the OPM modules using the same toolchain that
      was used to build DUNE. Otherwise, you can get strange ABI errors.


DEPENDENCIES FOR SUSE BASED DISTRIBUTIONS
-----------------------------------------

# repository containing prerequisites
sudo zypper ar http://download.opensuse.org/repositories/science/openSUSE_12.3/science.repo

# math libraries
sudo zypper in blas-devel lapack-devel suitesparse-devel superlu-devel

# utility libraries
sudo zypper in boost-devel tinyxml-devel

# tools necessary for building
sudo zypper in gcc gcc-c++ gcc-fortran cmake git doxygen

# DUNE libraries
sudo zypper in dune-common-devel dune-istl-devel

# Ensemble-based Reservoir Tools Eclipse utility module
git sudo zypper ar http://www.opm-project.org/packages/current/opensuse/12/opm.repo
sudo zypper in zlib-devel ert.ecl-devel

(to remove the repository, run `sudo zypper removerepo "Open Porous Media Initiative"`)

DEPENDENCIES FOR RHEL BASED DISTRIBUTIONS
-----------------------------------------

# packages necessary for building
sudo yum install make gcc-c++ gcc-gfortran cmake28 util-linux

# packages necessary for documentation
sudo yum install doxygen ghostscript texlive

# packages necessary for version control
sudo yum install git

# basic libraries necessary for both DUNE and OPM
sudo yum install boost-devel suitesparse-devel blas-devel lapack-devel

# libraries necessary for OPM
sudo yum install tinyxml-devel
sudo yum-config-manager --add-repo \
    http://www.opm-project.org/packages/current/redhat/6/opm.repo
sudo yum install libsuperlu3 ert.ecl-devel

# optional
sudo yum install dune-istl-devel


DEPENDENCIES FOR MACOS X
------------------------

You can build opm-core with Apple Xcode 4.6 or later, Ruby 1.9 or later
and the Homebrew port system:

# activate necessary repositories
brew tap homebrew/science
brew tap opm/opm

# libraries necessary for OPM
caffeinate brew install suite-sparse superlu ert.ecl
caffeinate brew install --with-c++11 boost tinyxml dune-istl


DOWNLOADING
-----------

For a read-only download:
git clone git://github.com/OPM/opm-core.git

If you want to contribute, fork OPM/opm-core on github.


BUILDING
--------

There are two ways to build the opm-core library.

1. As a stand-alone library.
In this setup we recommend creating an entirely separate directory
outside the directory containing the source code and doing the build
from that separate directory (termed "the build directory").  This
configuration is sometimes referred to as an "out-of-source build".

As an example, consider the following layout in which "opm-core" refers
to the directory containing the package source code as downloaded from
GitHub

    workspace
      |
      +-- build
      |
      +-- opm-core
      |     |
      |     +-- ...
      |     |
      |     +-- opm
      |     |
      |     +-- ...

We will configure a release-type (optimised) build using traditional
Unix Makefiles within the "build" directory.  The following command
configures the build

    cd path/to/build
    cmake ../opm-core -DCMAKE_BUILD_TYPE=Release

If you want to debug the library you should specify the build type
"Debug" instead of "Release" in the command above. This will disable
optimizations and make it easier to step through the code.

Building the software then amounts to typing

    make

in the top-level "build" directory; i.e., the directory from which we
invoked the "cmake" utility.  On a multi-core computer system you may
want to build the software in parallel (make(1)'s "job-server" mode) in
order to reduce the total amount of time needed to complete the build.
To do so, replace the above "make" command with

    make -j N

or, possibly,

    nice -20 make -j N

in which "N" is an integer that should typically not exceed the number
of cores in the system.

Once the library has been built, it can be installed in a central,
system-wide location (often in "/usr/local") through the command

    sudo make install


2. As a dune module.
 - Put the opm-core directory in the same directory
   as the other dune modules to be built (e.g. dune-commmon,
   dune-grid). Note that for Ubuntu you can install Dune
   from the ppa as outlined above.
 - Run dunecontrol as normal. For more information on
   the dune build system, see
   http://www.dune-project.org/doc/installation-notes.html


DOCUMENTATION
-------------

Efforts have been made to document the code with Doxygen.
In order to build the documentation, enter the command

 make doc

in the topmost directory.


REPORTING ISSUES
----------------

Issues can be reported in the Git issue tracker online at:

    http://github.com/OPM/opm-core/issues

To help diagnose build errors, please provide a link to a build log together
with the issue description.

You can capture such a log from the build using the `script' utility, e.g.:

    LOGFILE=$(date +%Y%m%d-%H%M-)build.log ;
	cmake -E cmake_echo_color --cyan --bold "Log file: $LOGFILE" ;
    script -q $LOGFILE -c 'cmake ../opm-core -DCMAKE_BUILD_TYPE=Debug' &&
    script -q $LOGFILE -a -c 'ionice nice make -j 4 -l 3' ||
    cat CMakeCache.txt CMakeFiles/CMake*.log >> $LOGFILE

The resulting file can be uploaded to for instance gist.github.com.
