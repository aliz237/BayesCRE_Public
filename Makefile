CXX = g++
CXXFLAGS = -std=c++14 -O2

main: fileio.o Relation.o Prob.o hash.o main.o
	$(CXX) -o $@ $^ $(CXXFLAGS) 
fileio.o: headers.h
Relation.o: headers.h 
Prob.o: headers.h
main.o: headers.h help.h
hash.o: headers.h

clean:
	rm -rf *.o *~ 
