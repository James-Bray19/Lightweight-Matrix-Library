# compiler settings
CC = gcc
CFLAGS = -Iinclude

# targets and dependencies
all: main

main: main.c lib/liblml.a
	$(CC) $(CFLAGS) $< -o $@ -Llib -llml

lib/liblml.a: src/lml.c
	$(CC) $(CFLAGS) -c $< -o $@

# clean up generated files
clean:
	rm -f main lib/liblml.a

