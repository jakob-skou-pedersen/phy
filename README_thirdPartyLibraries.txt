PHY LIBRARY THIRD PARTY LIBRARY DEPENDENCIES AND INSTALLATION
                       INSTRUCTIONS


OVERVIEW OF LIBRARY DEPENDENCIESOB

Essential:
* BOOST: Graph classes, various iterator classes, ublas
* BOOST numeric bindings: bindings to lapack, etc (using syev, etc)
* NTL: Number theoretic library (provides float type with extended exponent range)

Optional:
* opt++: Numeric optimization rutines
* newmat (v11): Linear algebra (Eigenvalue decomposition)
* bison++: grammar parser
* flex++: lexer 

Most functionality can used wihtout the optional libraries. 

The opt++ and newmat libraries are only needed for numerical
optimization routines. Most applications will not depend on these as
they currently only used for learning rate matrices describing the
substition process on phylogenetic trees.

Note that bison++ and flex++ are not needed to compile the code as
such. They are only needed to generate the code for the
newickParser. The FlexLexer header (FlexLexer.h), which is needed, is
provided with the code.

Below are given instructions for downloading source code and compiling
the third party libraries. Some libraries can probably be obtained as
packaged binaries.



INSTALLATION INSTRUCTIONS AND WEB-POINTERS

NECESSARY LIBRARIES:

=== BOOST ===

see http://www.boost.org/

 # fetch 
 wget http://sourceforge.net/projects/boost/files/boost/1.55.0/boost_1_55_0.tar.gz/download -O boost_1_55_0.tar.gz
 tar xvzf boost_1_55_0.tar.gz
 cd boost_1_55_0
 
 # compile and install
 ./bootstrap.sh
 ./bjam
 sudo ./bjam install
 

=== BOOST numeric bindings ===

see http://mathema.tician.de/node/391

 # fetch
 git clone http://git.tiker.net/trees/boost-numeric-bindings.git
 cd boost-numeric-bindings/

 # compile and install
 ./configure
 sudo make install
 # fix permissions
 sudo find /usr/local/include/boost-numeric-bindings/ -type d -exec chmod a+x {} \;
 sudo mv /usr/local/include/boost-numeric-bindings/boost/numeric/bindings /usr/local/include/boost/numeric/


=== NTL ===
see www.shoup.net/ntl/

 # fetch
 wget http://www.shoup.net/ntl/ntl-6.0.0.tar.gz
 tar xvzf ntl-6.0.0.tar.gz
 cd ntl-6.0.0/src

 # Compile and install
 ./configure
 make CFLAGS="-O2 -fPIC" # this takes 5 min...
 sudo make install



OPTIONAL LIBRARIES:

=== opt++ and newmat ===

====Overview====
opt++ comes bundles with newmat v11. We replace the bundled version
with the latest version and try to clean things up in a minimalistic
fashion.

====Homepages:====
* opt++: https://software.sandia.gov/opt++/
* newmt: http://www.robertnz.net/index.html

Registration is needed for download of opt++:
 https://software.sandia.gov/opt++/opt++_download.html

====compile and install ====
  # fetch
  tar xvzf tar xvzf optpp-2.4.tar.gz 
  cd optpp-2.4

  # compile and install
  configure --with-pic --includedir=/usr/local/include/optpp
  make
  sudo make install
  sudo chmod a+x /usr/local/include/optpp
  sudo cp include/OPT++_config.h /usr/local/include/optpp/

  # now the newmat include headers live in usr/local/include/optpp. We thus perform a
  # little hack to allow more standard calling of the include files
  # from code
  sudo ln -s /usr/local/include/optpp/ /usr/local/include/newmat

====Future improvements====

It would be nice to intall newmat indendently of opt++, which would
also allow newer versions to be used.

 
=== bison++ and flex++ === 
see http://code.google.com/p/flexpp-bisonpp/

 # fetch
 wget http://flexpp-bisonpp.googlecode.com/files/bisonpp-1.21-45.tar.gz
 tar xvzf bisonpp-1.21-45.tar.gz 
 cd bison++-1.21/

 # compile and install
 ./configure
 make
 sudo make install
 cd ..

 # fetch
 wget http://flexpp-bisonpp.googlecode.com/files/flexpp-2.3.8-45.tar.gz
 tar xvzf flexpp-2.3.8-45.tar.gz
 cd flex++-2.3.8-45/

 # compile and install
 ./configure 
 make
 sudo make install


For flex++, compilation problems have been reported in some cases,
which have been resolved by applying the fixes described here:

  http://code.google.com/p/flexpp-bisonpp/source/detail?r=26#

