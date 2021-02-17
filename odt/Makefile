all: code

code: code.o
	g++ -o rk4 code.o

code.o: code.cpp solver.h
	g++ -c code.cpp

clean:
	rm *.o
