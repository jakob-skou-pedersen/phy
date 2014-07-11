include ./config

all: NewickParserBase.cpp NewickParserBase.h NewickLexer.cpp
	$(MAKE) install

NewickParserBase.cpp NewickParserBase.h: newick.y
	bison++ -d -t -o NewickParserBase.cpp newick.y

NewickLexer.cpp: newick.l
	flex -+ -oNewickLexer.cpp  newick.l

install:
	cp -f NewickParserBase.h NewickParserBase.cpp NewickLexer.cpp ../

clean:
	rm -f *.o *~
	rm -f core*
	rm -f NewickLexer.cpp NewickParserBase.cpp 
	rm -f NewickParserBase.h


# # See http://code.google.com/p/flexpp-bisonpp/

# # The below urls also have packages and src
# # rpm bison++ package can be found at:
# # http://packages.debian.org/unstable/devel/bison++.html
# # http://ftp.de.debian.org/debian/pool/main/b/bison++/bison++_1.21.11.orig.tar.gz
# # http://www.mario-konrad.ch/file.php?id=16002
