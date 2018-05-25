all:
	gcc -Wall -g -O2 *.c -o main -lm

run:
	./main && gnuplot scen.plt