all:
	gcc -Wall -g -O0 *.c -o main -lm

run:
	./main && gnuplot "scen.plt"