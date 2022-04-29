CC = g++

.PHONY: clean

main: main.o
	$(CC) main.o -o main 

main.o: main.cpp matrix.hpp polynomial.hpp solver.hpp eigen.hpp approx.hpp
	$(CC) -c main.cpp

clean:
	rm *.o main
