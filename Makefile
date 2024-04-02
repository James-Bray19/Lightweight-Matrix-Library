CC = gcc
CFLAGS = -Wall -Wextra -Iinclude

.PHONY: all clean

all: example

example: examples/example.c src/matrixmagic.c
	$(CC) $(CFLAGS) $^ -o $@

clean:
	rm -f example
