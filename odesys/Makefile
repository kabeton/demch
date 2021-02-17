all: rk

rk: rk.o
	g++ -o rk rk.o

rk.o: rk.cpp rk.h
	g++ -c rk.cpp

clean:
	rm *.o
