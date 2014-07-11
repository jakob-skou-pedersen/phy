INSTALLATION INSTRUCTIONS FOR THE PHY LIBRARY

The phy library provides data structures, algorithms, and io routines
for discrete factor graphs (DFGs), phylogenetic substitution models,
and stochastic context-free grammars (SCFGs). A number of applications
based on this functionality (EvoFold, dfgEval, and dfgTrain, etc) are
included in the src directory. 

The library makes use of namespaces. Nearly all the functionality is
included in the phy namespace. The grammar namespace is nested within
phy and provides the SCFG related functionality.

The main documentation of the library is currently provided by
comments within the source code. Header files often start with a short
description of the organization and functionality of included classes,
etc.. In addition, each class and function is preceeded by a short
description.  

The phy library relies on third party libraries for different types of
functionality, including: iterator classes, unnit testing, c++
bindings to numerical functions, etc (BOOST); matrix algebra (ublas);
floating point types with extended exponent ranges (NTL); numerical
optimization (Newmat and OPT++); parser generation (Flex and Bison).

Not all of these libraries are needed for all functionality. The OPT++
and Newmat libraries can optionally be left out using a preprocessor
directive (NO_OPTPP) and Flex and Bison can be left out using a
similar directive (NO_FLEX). In this case the numerical optimization
routines are not provided and phylogenetic trees in the Newick format
cannot be parsed. This functionality is currently only used for
phylogenetic analysis.

The installation of the third party libraries is descibed in
README_thirdPartyLibraries.txt.

Below are given the compilation and installation instructions under
Linux and Mac. 


LINUX:

Fist install the third party libraries as described in
README_thirdPartyLibraries.txt

If Autoconf, Automake, and Libtool are not present, these need to be
installed first. If using apt-get for package management, this can be
done by:

  sudo apt-get install libtool automake autoconf

To generate configuration files and makefiles, run: 

  autoreconf --force --install

The flags used for compilation are set at the configure step. We thus
point to the location(s) of header files and library files at this
stage. In the below we assume these are found under /usr/local/. Note
the use of -isystem instead of -I, which avoids warnings from third
party libraries. By default the library and applications are installed
under /usr/local/. This can be changed using the --prefix flag.

As the configuation depends on library use, we provide two different
versions:

A) For production (fastest running applications):

  ./configure CXXFLAGS="-Wall -O2 -DNDEBUG -isystem /usr/local/include -D NO_OPTPP -D NO_FLEX"

B) For development (including debugging information):

  ./configure CXXFLAGS=" -Wall -ggdb3 -O0 -isystem /usr/local/include -D NO_OPTPP -D NO_FLEX" 

To include OPT++ and Flex / Bison libraries leave out "-D NO_OPTPP -D
NO_FLEX" from the above.


Compilation, checking (running a unit test suite), and installation of
phy library files:

  make
  make check
  make install # as root

Compilation and installation of executables 

  cd ./phy/src 
  make 
  make install # as root


MAC: 

Installation and compilation on MAC should generally follow the Linux
instructions. You may need to point the linker to the location of the
shared libraries if they are in non-standard location, such as
/opt/local/lib:

  export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:/opt/local/lib

