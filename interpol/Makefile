all: code

code: code.o 
	g++ -o inter code.o 

code.o: code.cpp poly.h spline.h
	g++ -c code.cpp

clean:
	rm *.o
