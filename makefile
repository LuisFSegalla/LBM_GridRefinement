all:
	gcc -o main2D main2D.c q2d9_library.c -lm

roda:
	./main2D
	
limpa:
	rm -f *.o

limpaTudo:
	rm -f *.0 *.txt
