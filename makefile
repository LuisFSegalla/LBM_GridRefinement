all:
	gcc -g -o main2D.out main2D.c q2d9_library.c -lm

roda:
	./main2D
	
limpa:
	rm -f *.o

limpaTudo:
	rm -f *.0 *.txt results/*.txt
