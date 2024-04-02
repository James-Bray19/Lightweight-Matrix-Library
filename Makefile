# compiler settings
CC = gcc
CFLAGS = -Iinclude

# targets and dependencies
all: example

example: main.c lib/libmatrixmagic.a
	$(CC) $(CFLAGS) $< -o $@ -Llib -lmatrixmagic

lib/libmatrixmagic.a: src/matrixmagic.c
	$(CC) $(CFLAGS) -c $< -o $@

# clean up generated files
clean:
	rm -f example lib/libmatrixmagic.a

