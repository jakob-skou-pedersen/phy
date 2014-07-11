calc:
	bison++ -d -o parser.cc exp.y 
	flex++ -I -olexer.cc exp.l
	g++ -o calc -Wall *.cc -lfl -s 


clean:
	rm -f *.o parser.h parser.cc lexer.cc calc
