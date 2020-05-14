
# Basic parameters
CC     = gcc 
CFLAGS = -std=c99 -Wall -Wextra -pedantic -O3
OBJS   = MedianSort.o

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

test: $(OBJS)
	$(CC) $(CFLAGS) $^ -o $@
	./test > log.txt
	gnuplot < format.plt
	/bin/rm -rf *.o *~
	/bin/rm -rf test
